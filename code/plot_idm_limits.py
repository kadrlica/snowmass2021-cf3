#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
__editor__ = "Ethan O. Nadler, Sidney Mau"
import copy
import yaml
import numpy as np

import pylab as plt
import matplotlib
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import interp1d

import plot
from plot import plot_limit, plot_limit_fill, plot_limit_patch, get_mass_limit, plot_text
from plot import custom_blues

rcParams = matplotlib.rcParams
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
rcParams["xtick.direction"] = 'in'
rcParams["ytick.direction"] = 'in'
rcParams["xtick.major.size"] = 6
rcParams["ytick.major.size"] = 6
rcParams["xtick.minor.size"] = 3
rcParams["ytick.minor.size"] = 3
rcParams["xtick.major.width"] = 1.0
rcParams["ytick.major.width"] = 1.0
rcParams["xtick.minor.width"] = 0.8
rcParams["ytick.minor.width"] = 0.8
rcParams["xtick.major.pad"] = 5
rcParams["ytick.major.pad"] = 5

matplotlib.rcParams["legend.handlelength"] = 2.2

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--arrows',action='store_true')
parser.add_argument('--summary',action='store_true')
parser.add_argument('--project',action='store_true')
args = parser.parse_args()

draw_arrows = args.arrows
draw_summary = args.summary
draw_project = args.project

filename = plot.get_datafile('idm_limits.yaml')
limits = yaml.load(open(filename))

f,ax = plt.subplots(figsize=(12,8))
ax.set_xscale('log'); ax.set_yscale('log')

# Direct Detection
plot_limit_patch(limits['xenon1t_cresst3_dd'],alpha=1.0,color='0.8',zorder=2)
#plot_limit(limits['xenon1t_cresst3_dd'],alpha=1.0,color='0.8',zorder=0)
plot_limit_fill(limits['mahdawi2019_xqc_e002'],alpha=0.7,color='0.8',zorder=2);
plot_limit(limits['mahdawi2019_xqc_e002'],color='k',lw=1.5,ls='--',zorder=10);
plot_text(limits['mahdawi2019_xqc_e002'])
ax.annotate('Direct Detection',(1e0,1e-35),color='k',fontsize=22)

if draw_arrows:
    arrow = FancyArrowPatch((2e0, 1e-38), (2e-1, 1e-39), mutation_scale=50, fc='0.8', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((3e1, 1e-45), (4.5e1, 3e-49), mutation_scale=50, fc='0.8', ec='none')
    ax.add_patch(arrow)

# Neutrino Fog
#plot_limit_patch(limits['neutrino_fog_n2'],zorder=999,edgecolor='none')
#plot_limit_patch(limits['neutrino_fog_n2_ext'],zorder=999,edgecolor='blue',alpha=0.5)
#plot_limit_patch(limits['ohara2021_neutrino_floor'],zorder=1)
plot_limit_patch(limits['nuetrino_xe_si_2'],zorder=1)

ax.annotate('Neutrino Fog',(2e-1,1e-47),color='k',fontsize=22)

# Cosmic Probes
if draw_summary:
    plot_limit_patch(limits['cosmic_idm'],color='g',zorder=10)
    plot_limit_patch(limits['cosmic_idm'],color='w',alpha=1.0,zorder=1)
    plot_limit(limits['cosmic_idm'],color='g',alpha=1.0,zorder=10)
    ax.annotate("CMB and Small-Scale Structure",(1e-4,1e-28),color='k',fontsize=16)
    ax.annotate("Small-Scale Structure",(3e-6,1e-36),color='k',fontsize=16,rotation=90)
    #ax.annotate('Cosmic Probes',(1e-3,1e-27),color='k',fontsize=22)
    #ax.annotate('Cosmic Probes',(2e-5,1e-31),color='g',fontsize=30)
    ax.annotate(r'{\bf Cosmic Probes}',(2e-6,1e-26),color='k',fontsize=30,zorder=999)
else:
    # IDM
    plot_limit_fill(limits['rogers2022_lya'],alpha=0.7,zorder=10); plot_text(limits['rogers2022_lya'])
    plot_limit_fill(limits['nadler2021_dsphs'],alpha=0.7,zorder=10); plot_text(limits['nadler2021_dsphs'])
    plot_limit_fill(limits['gluscevic2018_planck'],alpha=0.7,zorder=10); plot_text(limits['gluscevic2018_planck'])
    # WDM
    plot_limit_fill(limits['nadler2022_dsphs_sl'],zorder=10); plot_text(limits['nadler2022_dsphs_sl'])
    plot_limit_fill(limits['gilman2019_sl'],zorder=10); plot_text(limits['gilman2019_sl'])
    plot_limit_fill(limits['irsic2017_lya'],zorder=10); plot_text(limits['irsic2017_lya'])
    #ax.annotate('Cosmic Probes',(2e-5,1e-31),color='g',fontsize=30)
    ax.annotate(r'{\bf Cosmic Probes}',(2e-6,1e-26),color='k',fontsize=30,zorder=999)

#mass = np.logspace(-6.0,-2.5,100)
#sigma[mass > 9.7e-6] = 3.67244e-30
#sigma[mass > 8.3e-5] = 6.27456e-30
#sigma[mass > 3.7e-4] = 8.61084e-30
#sigma[mass > 1.2e-3] = 1.12554e-29
#sigma[mass > 3.4e-3] = 1.40128e-29
#sigma[mass > 8.0e-3] = 1.66167e-29
mass,sigma = get_mass_limit(limits['cosmic_idm'])
interp = interp1d(mass[sigma > 1e-46],sigma[sigma > 1e-46])

mass = np.logspace(np.log10(9.7e-6),-2.5,100)
sigma = np.logspace(-49,-30,100)

for x in range (75, 100):
    sigma = [sigma[0],interp(mass[x])]
    #ax.fill_betweenx(sigma, mass.min(), mass[x], facecolor='green', alpha=0.02, zorder=0)
    ax.fill_between(mass[:x], sigma[0], sigma[1], facecolor='green', alpha=0.02, zorder=0)
ax.annotate(r'Big Bang Nucleosynthesis',(1e-3,1e-36),fontsize=14,rotation=90)


if draw_arrows:
    arrow = FancyArrowPatch((9.6e-6, 1e-37), (5e-5, 1e-37), mutation_scale=70, fc='g', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((1e-2, 2.6e-29), (1.3e-2, 2e-32), mutation_scale=70, fc='g', ec='none')
    ax.add_patch(arrow)
elif draw_project:
    plot_limit(limits['cosmic_project'],zorder=10)

#Ticks, axes
ymin,ymax=1e-50,2e-24
xmin,xmax = 1e-6,1e2

ax.set_xlim(xmin=xmin,xmax=xmax)
ax.set_ylim(ymax=ymax,ymin=ymin)

#ax.set_xticks([1e-6,1e-4,1e-2,1e0,1e2])
#ax.set_xticklabels([r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$',r'$10$',r'$10^{2}$'],fontsize=22)

xticks = 10.**np.arange(-6,3)
xticklabels = [r'$10^{%g}$'%i for i in np.arange(-6,3)]

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels,fontsize=22)


ax.set_yticks([1e-45,1e-40,1e-35,1e-30,1e-25])
ax.set_yticklabels([r'$10^{-45}$',r'$10^{-40}$',r'$10^{-35}$',r'$10^{-30}$',r'$10^{-25}$'],fontsize=22)

ax.set_xlabel('Dark matter mass (GeV/$c^2$)',fontsize=26,labelpad=8)
ax.set_ylabel('Dark matter-nucleon cross section (cm$^2$)',fontsize=26,labelpad=12)

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)


plt.tight_layout()

outfile = 'idm_limits.pdf'
if draw_summary:
    plt.savefig(outfile.replace('.pdf','_summary.pdf'))
    plt.savefig(outfile.replace('.pdf','_summary.png'))
else:
    plt.savefig(outfile)
    plt.savefig(outfile.replace('.pdf','.png'))
