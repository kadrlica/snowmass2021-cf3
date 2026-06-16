#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import copy
import pylab as plt
import matplotlib
import numpy as np
import yaml

import plot
from plot import plot_limit, plot_limit_fill, plot_limit_patch, get_mass_limit
from plot import custom_blues

matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['legend.handlelength'] = 2.2

def m22_nu(m22):
    return 1e-6 * (m22/40)

def nu_m22(nu):
    return  40 * nu / 1e-6

filename = plot.get_datafile('fdm_limits.yaml')
limits = yaml.load(open(filename))

xticks=nu_m22(np.array([1e-6,1e-3,1,1e3,1e6]))*1e-22
xticklabels=[r'$\mu$Hz',r'mHz',r'Hz',r'kHz',r'MHz']
#xticklabels=[r'$10^{-6}$',r'$10^{-3}$',r'$1$',r'$10^{3}$',r'$10^{6}$']

fig,ax = plt.subplots(figsize=(8,6))
ax.set_yscale('log')
ax.set_xscale('log')
ax2 = ax.twiny()

plot_limit(limits['coleman2018_magis100'])
ax.annotate('MAGIS-100',xy=(2e-18,3e-28),xycoords='data',fontsize=18)


plot_limit_fill(limits['graham2016_EP'])
plot.plot_text(limits['graham2016_EP'])
ax.annotate('Torsion Balance Tests \n of the Equivalence Principle',
            xy=(1e-17,5e-23),xycoords='data',fontsize=18)


plot_limit_fill(limits['schutz2020_streams'])
schutz = get_mass_limit(limits['schutz2020_streams'])
ax.loglog(schutz[0],schutz[1],c=custom_blues[4],lw=3,zorder=999,label=r'Current Sensitivity')

delve = get_mass_limit(limits['delve'])
ax.loglog(delve[0],delve[1],c='r',ls='--',lw=3,zorder=999,
          label=r'DECam Sensitivity')

future = get_mass_limit(limits['lsst_streams'])
ax.fill_between(future[0],future[1],1,color='r',alpha=0.3)
ax.loglog(future[0],future[1],c='r',lw=3,zorder=999,
          label=r'LSST Sensitivity')

#ax.annotate('Stellar \n Streams',xy=(1e-20,1.5e-26),xycoords='data',
#            rotation=90,fontsize=20,fontweight='bold')

ax2.set_xscale('log')
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticklabels)

ax2.set_xlim(1e-22,1e-11)
ax.set_xlim(1e-22,1e-11)
ax.set_ylim(1e-32,1e-20)
ax.set_xlabel(r'Dark Matter Mass (eV)',fontsize=18)
ax.set_ylabel(r'Dark Matter Vector Coupling ($g_{B-L}$)',fontsize=18)
ax.legend(loc='lower right',prop={'size':16},frameon=False)

plt.subplots_adjust(top=0.90,bottom=0.12,left=0.15)

plt.savefig('fdm_limits_gBL.pdf')
plt.savefig('fdm_limits_gBL.png')
plt.ion()

