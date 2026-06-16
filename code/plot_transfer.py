#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Ethan O. Nadler"
__editor__ = "Sidney Mau"

# Imports

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import plot
from plot import custom_blues

matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['legend.handlelength'] = 2.2

custom_cet = {}
custom_cet['wdm_lovell_new'] = ["#FFBB99", "#FF9966", "#FF7733", "#FF5500", "#CC4400", "#993300", "#662200","#000000"]
custom_cet['wdm_lovell'] = ["#c3ff14","#33dc14","#2ec611"]
custom_cet['wdm_lovell_alt'] = ['brown']
custom_cet['wdm_schneider'] = ['gold']
custom_cet['idm'] = ["#99DDFF","#66CCFF", "#33BBFF", "#00AAFF", "#0088CC", "#006699", "#004466", "#000000"]
custom_cet['fdm_du'] = ["#fc5eff","#e330ff","#b11eff"]
custom_cet['fdm_schive'] = ['brown']

custom_blues = ["#99DDFF","#66CCFF", "#33BBFF", "#00AAFF", "#0088CC", "#006699", "#004466", "#000000"]
custom_blues_complement = ["#FFBB99", "#FF9966", "#FF7733", "#FF5500", "#CC4400", "#993300", "#662200","#000000"]

# Cosmological Parameters

h = 0.7
omega_m = 0.286
rho_crit = 1.36*10**11*omega_m #msun/Mpc^3

# Input Parameters

Mmin = {}
Mmin['DES'] = 3.2e8
Mmin['SDSS'] = 5.4e8

illustrative_params = {}
illustrative_params['wdm'] = 5e8
illustrative_params['fdm'] = 10

ruled_out_params = {}
ruled_out_params['wdm_lovell_new'] = 9.4e7
ruled_out_params['wdm_lovell'] = 3.8e7
ruled_out_params['fdm_du'] = 28
ruled_out_params['fdm_schive'] = 91

# Transfer functions for T^2(k)

def khm(mwdm):
    lambda_fs = (0.049*(mwdm**(-1.11))*((omega_m/0.25)**(0.11))*((h/0.7)**1.22))
    lambda_hm = 13.93*lambda_fs
    k_hm = 2*np.pi/lambda_hm
    return k_hm

def mwdm(Mhm):
    lambda_hm = 2.*(((3./4)*Mhm/(np.pi*rho_crit))**(1./3))
    lambda_fs = lambda_hm/13.93
    mwdm = (h*lambda_fs/(0.049*((omega_m/0.25)**(0.11))*((h/0.7)**1.22)))**(1./-1.11)
    return mwdm

def khm2(Mhm):
    return khm(mwdm(Mhm))

def T2_wdm(k,mwdm):
    nu = 1.12
    lambda_fs = (0.049*(mwdm**(-1.11))*((omega_m/0.25)**(0.11))*((h/0.7)**1.22))
    alpha = lambda_fs
    transfer = (1+(alpha*k)**(2*nu))**(-10./nu)
    return transfer

#Load IDM transfer function

filename0 = plot.get_datafile('lcdm_z1_pk.dat')
data0 = np.loadtxt(filename0)
ks = data0[:,0]
pk0 = data0[:,1]

filename = plot.get_datafile('n0_m0.001_s4.83293023857e-29_z1_pk.dat')
data = np.loadtxt(filename)
pk = data[:,1]

def kj(mfdm):
    return 9.*((mfdm)**0.5)/h

def xj(k,mfdm):
    return 1.61*((mfdm)**(1./18))*k/kj(mfdm)

def T2_fdm(k,mfdm):
    return (np.cos(xj(k,mfdm)**3)/(1+xj(k,mfdm)**8))**2

# SHMFs for f_{DM}(M, \xi_{DM})

def cdm(Mpeak):
    return 3.26*10**-5 * (1./(2.57*10**7))**(-1.9) * (Mpeak)**(-1.9+1)

mpeak = np.logspace(7,10,100)
avg = cdm(mpeak)

def wdm_lovell_new(Mpeak,Mhm):
    return (1+(4.2*Mhm/Mpeak)**2.5)**-0.2

def wdm_lovell(Mpeak,Mhm):
    return (1+2.7*Mhm/Mpeak)**-0.99

def wdm_lovell_alt(Mpeak,Mhm):
    return (1+Mhm/Mpeak)**-1.3

def wdm_schneider(Mpeak,Mhm):
    return (1+Mhm/Mpeak)**-1.16

def fdm_du(avg,mpeak,m22):
    M1 = 4.7*(m22**-1.5)*10**8
    M2 = 2.0*(m22**-1.6)*10**8
    f1 = 0.014*(m22**1.5)*np.exp(-np.log(mpeak/M1)**2/1.4)
    f2 = (1. + (mpeak/M2)**-0.72)**(-10./0.72)
    return f1 + f2*avg

def fdm_schive(Mpeak,m22,alpha=-1.1,beta=-2.2):
    M0 = 1.6*10**10*(m22**(-4./3.))
    return (1+(Mpeak/M0)**alpha)**beta

# Plots

fig = plt.figure(figsize=(16,7))

#Transfer function panel
ax = fig.add_subplot(121)
ax2 = ax.twiny()
ax.set_xscale('log')

ks_array = np.linspace(0.1,250,1000)

ax.semilogx(ks_array,T2_wdm(ks_array,mwdm(ruled_out_params['wdm_lovell'])),
                  lw=3,linestyle='-',c=custom_cet['wdm_lovell_new'][3],
                  label=r'Warm Dark Matter')

ax.semilogx(3.4*ks,pk/pk0,
            alpha=1.0,c=custom_cet['idm'][3],lw=3,ls='--',
            label=r'Interacting Dark Matter')

ax.semilogx(ks_array,T2_fdm(ks_array,ruled_out_params['fdm_du']),
            linestyle='-',lw=3,c=custom_cet['fdm_du'][2],
            label=r'Fuzzy Dark Matter')

ax.semilogx(ks_array[ks_array<150],ks_array[ks_array<150]/ks_array[ks_array<150],'k--',lw=1,zorder=1)

ax.legend(loc=3,fontsize=20,frameon=False)

ax.set_xlim(2,150)
ax.set_ylim(0.0005,1.05)
ax.set_xlabel(r'$k$ [$h\ \rm{Mpc}^{-1}$]',fontsize=28)
ax.set_ylabel(r'$T^2(k)$',fontsize=28,labelpad=12)
ax.set_xticks([2,5,10,30,50,100])
ax.set_xticklabels([r'$2$',r'$5$',r'$10$',r'$30$',r'$50$',r'$100$'],fontsize=22)
ax.set_yticks([0.25,0.5,0.75,1.0])
ax.set_yticklabels([0.25,0.5,0.75,1.0],fontsize=22)

ax2.set_xscale('log')
ax2.set_xlim(2,150)
ax2.set_xticks([khm2(10**12),khm2(10**11),khm2(10**10),khm2(10**9),khm2(10**8),khm2(10**7)])
ax2.set_xticklabels([r'10$^{12}$',r'10$^{11}$',r'10$^{10}$',r'10$^{9}$',r'10$^{8}$',r'10$^{7}$'],fontsize=22)
ax2.set_xlabel(r'$M$ [$M_{\rm{\odot}}$]',fontsize=28,labelpad=12)

#SHMF panel
ax = fig.add_subplot(122)

ax.loglog(mpeak,wdm_lovell(mpeak,ruled_out_params['wdm_lovell']),
           c=custom_cet['wdm_lovell_new'][3],lw=3,
           label=r'WDM $(m_{\mathrm{WDM}}=5.0\ \mathrm{keV})$')

ax.loglog(0.95*mpeak,1.05*wdm_lovell(mpeak,ruled_out_params['wdm_lovell']),
           c=custom_cet['idm'][3],lw=3,ls='--',alpha=1.0,
           label=r'IDM $(\sigma_{0}=3\times 10^{-30}\ \mathrm{cm}^2)$')

ax.loglog(mpeak,fdm_du(avg,mpeak,ruled_out_params['fdm_du'])/avg,
           c=custom_cet['fdm_du'][2],lw=3,
           label=r'FDM $(m_{\phi}=2.8\times 10^{-21}\ \mathrm{eV})$')

ax.loglog(Mmin['DES']*np.ones(10),np.linspace(0.85e-1,1.25,10),c='gray',ls='--',zorder=1)
ax.text(2.15e8,0.4,r'Minimum Halo Mass',color='k',alpha=0.75,fontsize=22,rotation=90)

ax.set_xlim(6e6,1.25e10)
ax.set_ylim(0.85e-1,1.25)
ax.set_xlabel(r'$M\ [M_{\mathrm{\odot}}]$',fontsize=28,labelpad=8)
ax.set_ylabel(r'$f_{\mathrm{DM}}(M,\xi_{\mathrm{DM}})$',fontsize=28,labelpad=12)
ax.set_xticks([1e7,1e8,1e9,1e10])
ax.set_xticklabels([r'$10^7$',r'$10^8$',r'$10^9$',r'$10^{10}$'],fontsize=22)
ax.set_yticks([0.1,0.25,0.5,0.75,1.0])
ax.set_yticklabels([r'$0.1$',r'$0.25$',r'$0.5$',r'$0.75$',r'$1.0$'],fontsize=22)

# ax2 = ax.twiny()
# ax2.set_xlim(10**7,10**11)
# ax2.set_xscale('log')
# ax2.set_xticks([26242185.4338,197696964.011,1489361077.71,11220184543,84527884516])
# ax2.set_xticks([107645701.875,810954884.192,6109373739.43,46025306975.3],minor=True)
# ax2.set_xticklabels([r'$10^7$',r'$10^8$',r'$10^9$',r'$10^{10}$',r'$10^{11}$'],fontsize=18)
# ax2.set_xticklabels([r'',r'',r'',r''],minor=True)
# ax2.set_xlabel(r'$\mathcal{M}_{\rm{vir}}$ [$M_{\rm{\odot}}$]',fontsize=22,labelpad=8)

plt.tight_layout()
plt.gcf().subplots_adjust(wspace=0.25)
plt.savefig('transfer.pdf')

#plt.show()
