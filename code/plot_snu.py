#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
__editor__ = "Ethan O. Nadler, Sidney Mau"
import copy
import pylab as plt
import numpy as np
import yaml
import matplotlib

from matplotlib.path import Path
from matplotlib.patches import PathPatch

import plot
from plot import plot_limit, plot_limit_fill, plot_limit_patch, get_mass_limit
from plot import custom_blues

matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['legend.handlelength'] = 2.2

h = 0.6727
omega_m = (0.1199)/(h**2)

def mwdm2msnu(mwdm):
    # Bullock & Boylan-Kolchin
    # https://arxiv.org/pdf/1707.04256.pdf
    msnu = 3.9 * (mwdm)**1.294 * (omega_m * h**2 / 0.1225)**(-1/3.)
    return msnu

def logline(x,m,b):
    return 10**(m*np.log10(x) + b)

def overprod(x):
    # Boyarsky 2009
    #https://arxiv.org/pdf/0901.0011.pdf
    # Adjusted to Schneider 2016
    m = -1.7847
    #b = -7.3169 # original
    b = -7.192 # adjusted
    return logline(x,m,b)

def bbn(x):
    # Boyarsky 2009
    #https://arxiv.org/pdf/0901.0011.pdf
    m = -1.4257
    b = -11.9965
    return logline(x,m,b)

def nuMSSM(x):
    # Boyarsky 2009
    #https://arxiv.org/pdf/0901.0011.pdf
    m = -1.40098
    b = -11.51094
    return logline(x,m,b)

filename = plot.get_datafile('snu_limits.yaml')
limits = yaml.load(open(filename))

fig,ax = plt.subplots(figsize=(8,6))
ax.set_yscale('log')
ax.set_xscale('log')

plot_limit_fill(limits['tremaine_gunn'])
plot_limit(limits['tremaine_gunn'])

#plot_limit_fill(limits['merle2017_katrin'])
#plot_limit(limits['merle2017_katrin'])

#plot_limit_fill(limits['merle2017_xray'])
#plot_limit(limits['merle2017_xray'])

#plot_limit_fill(limits['merle2017_combined_xray'])
#plot_limit(limits['merle2017_combined_xray'])

#plot_limit_fill(limits['cherry2017_xray'])
#plot_limit(limits['cherry2017_xray'])

plot_limit_fill(limits['combined_xray'])
plot_limit(limits['combined_xray'])

#plot_limit(limits['horiuchi2013_m31'])
#plot_limit(limits['perez2016_nustar_gc'])

#plot_limit(limits['ng2019_nustar_m31'])
#plot_limit(limits['roach2019_nustar_bulge'])

#plot_limit(limits['boyarsky2007_integral'])
#plot_limit(limits['boyarsky2006_heao1_mw'])

#plot_limit(limits['cherry2017_maximal'])
#plot_limit(limits['cherry2017_sdss_m31'])
#plot_limit(limits['cherry2017_sdss_allsky'])
#plot_limit(limits['cherry2017_sdss_des'])

#plot_limit_fill(limits['dessert2020_xmm_lines'])

#plot_limit_fill(limits['schneider2016_sats'])
current = get_mass_limit(limits['schneider2016_sats'])

ax.text(45,5e-9,r'DM Overproduction',rotation=-22,fontsize=22,ha='center',  va='top')
ax.text(2.75,2.5e-12,r'DM Underproduction',rotation=-17,fontsize=22,ha='center',  va='top')

#Bounds
m = np.array([np.min(current[0]),300])
ax.loglog(m,nuMSSM(m),ls='-',c='black',lw=2.5)

#plot_limit_fill(limits['merle2017_overprod'])
#plot_limit(limits['merle2017_overprod'])
#m,l = get_mass_limit(limits['merle2017_overprod'])

ax.fill_between(m,np.ones(len(m))*1e-15,nuMSSM(m),facecolor='black',alpha=0.1)

plot_limit_fill(limits['schneider2016_overprod'])

m = np.array([np.min(current[0]),300])
ax.loglog(m,overprod(m)*1.0,lw=2.5,c='k',ls='-',zorder=1)

#plot_limit(limits['boyarsky2009_bbn'])

#SDSS
msnu0 = 10
mix0 = overprod(msnu0)

sdss = np.copy(current)
msnu = mwdm2msnu(2.5)
mix = overprod(msnu)
sdss[0] *= msnu/msnu0
sdss[1] *= mix/mix0
ind = np.logical_and(sdss[1]>nuMSSM(sdss[0]),sdss[1]<overprod(sdss[0]))

ax.loglog(np.array([6.,6.]),
          np.array([nuMSSM(sdss[0])[np.argmin(np.abs(sdss[0]-6.))],overprod(sdss[0])[np.argmin(np.abs(sdss[0]-6.))]]),
          ls='-',c=custom_blues[3],lw=3,zorder=0,label=r'SDSS Dwarfs')

ax.fill_between(np.linspace(0.5,6.0,100),nuMSSM(np.linspace(0.5,6.0,100)),overprod(np.linspace(0.5,6.0,100)),
                edgecolor="None",
                facecolor=custom_blues[3],
                alpha=0.35)

#DES
msnu0 = 10
mix0 = overprod(msnu0)

decam = np.copy(current)
msnu = mwdm2msnu(6)
mix = overprod(msnu)
decam[0] *= msnu/msnu0
decam[1] *= mix/mix0
ind2 = np.logical_and(decam[1]>nuMSSM(decam[0]),decam[1]<overprod(decam[0]))

ax.loglog(np.array([25.0,25.0]),
          np.array([nuMSSM(decam[0])[np.argmin(np.abs(decam[0]-25.0))],overprod(decam[0])[np.argmin(np.abs(decam[0]-25.0))]]),
          ls='-',c='r',lw=3.5,zorder=0,label=r'DES Dwarfs')

ax.fill_between(np.linspace(6.0,25.,100),nuMSSM(np.linspace(6.0,25.,100)),overprod(np.linspace(6.0,25.,100)),
                edgecolor="None",
                facecolor='r',
                alpha=0.35)

#lsst = np.copy(current)
#msnu = mwdm2msnu(20)
#mix = overprod(msnu)
#lsst[0] *= msnu/msnu0
#lsst[1] *= mix/mix0
#ax.loglog(lsst[0],lsst[1],ls='-',c='r',lw=3,zorder=0,label=r'LSST Sensitivity')
#
#x = np.hstack([current[0][1:-1],lsst[0][1:-1][::-1]])
#y = np.hstack([current[1][1:-1],lsst[1][1:-1][::-1]])
#path = Path(np.vstack([x,y]).T)
#patch = PathPatch(path,facecolor='r',edgecolor='none',alpha=0.3)
#ax.add_artist(patch)
##ax.fill_betweenx(lsst[1][1:],current[0][:-1],lsst[0][1:],color='r',alpha=0.3,hatch='/')
#
##lsst = np.copy(current)
##xfactor = 4
##lsst[0] *= 16
##lsst[1] /= overprod(lsst[0])
##ax.loglog(lsst[0],lsst[1],ls='-',c='r',lw=3,zorder=0,label=r'LSST Sensitivity')
#
#
##future = get_mass_limit(limits['lsst_streams'])
##plot.plot_limit_patch(limits['lsst_streams'],zorder=1)
##ax.loglog(np.nan,np.nan,c='r',lw=3,zorder=999,label=r'LSST Sensitivity')
##ax.fill_between(current[1:-1],current[1:-1],1,color='r',alpha=0.3,hatch='/')
###ax.fill_between(future[0],future[1],1,color='r',alpha=0.3,hatch='/')
##ax.loglog(future[0],future[1],c='r',lw=3,zorder=999,label=r'Future Sensitivity')
#
#
##m = np.array([0.5,100])
##ax.loglog(m,nuMSSM(m),ls='--',label='nuMSSM')
##ax.loglog(m,bbn(m),ls='--',label='BBN')
##ax.loglog(m,overprod(m),ls='--',label='Overprod')
#
##plot_limit_fill(limits['merle2017_overprod'])
##plot_limit(limits['merle2017_overprod'])
##m,l = get_mass_limit(limits['merle2017_overprod'])
#
#plot_limit_fill(limits['schneider2016_overprod'])
#plot_limit(limits['schneider2016_overprod'])
#m,l = get_mass_limit(limits['schneider2016_overprod'])
#ax.loglog(m,l,lw=3.0,ls='-',color='0.2')
#
##plot_limit(limits['boyarsky2009_bbn'])

# Boyarsky 2014
x0, y0 = 7, 5.11656e-11
ylo,yhi =2.17501e-11,1.98250e-10
yerr = np.array([[y0-ylo, yhi-y0]]).T
ax.errorbar(x0,y0,yerr=yerr,fmt='s',color='k',ms=5)

# Ng et al. 2019 (scraping the point to double check)
#https://arxiv.org/abs/1901.01262
x0,y0 = 7.06, 5.00841e-11
ylo,yhi=2.03849e-11,1.95565e-10
yerr = np.array([[y0-ylo, yhi-y0]]).T
#ax.errorbar(x0,y0,yerr=yerr,fmt='s',color='r',markersize=3)


#ax.set_xlim(0.3,50)
#ax.set_ylim(5e-13,1e-2)

#Ticks, axes
ax.set_xlim(0.3,250)
ax.set_ylim(1e-15,1e-6)
#ax.set_yticks([1e-12,1e-9,1e-6,1e-3])

ax.set_xticks([0.5,1,5,10,50,100])
ax.set_xticklabels([r'$0.5$',r'$1$',r'$5$',r'$10$',r'$50$',r'$100$'],fontsize=22)

#ax.set_yticks([1e-15,1e-13,1e-11,1e-9,1e-7,1e-5,1e-3])
ax.set_yticks([1e-15,1e-13,1e-11,1e-9,1e-7])
ax.set_yticklabels([r'$10^{-15}$',r'$10^{-13}$',r'$10^{-11}$',r'$10^{-9}$',r'$10^{-7}$',r'$10^{-5}$',r'$10^{-3}$'],fontsize=22)

ax.set_xlabel(r'$m_s$ (keV)',fontsize=26,labelpad=8)
ax.set_ylabel(r'$\sin^2(2\theta)$',fontsize=26,labelpad=12)
ax.set_title(r'Sterile Neutrino WDM',fontsize=28,y=1.01)
ax.legend(prop={'size':18},frameon=True,loc=1,handletextpad=0.5)

plt.subplots_adjust(left=0.15,top=0.95,bottom=0.12,right=0.95)

plt.tight_layout()

plt.savefig('snu_limits.pdf')
plt.savefig('snu_limits.png')
plt.ion()
