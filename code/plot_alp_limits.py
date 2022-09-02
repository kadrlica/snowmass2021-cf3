#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
__editor__ = "Ethan O. Nadler, Sidney Mau"
import os
import copy
import yaml
import numpy as np
from scipy.interpolate import interp1d
import pylab as plt
import matplotlib
from matplotlib.patches import FancyArrowPatch


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
rcParams["xtick.major.pad"] = 7
rcParams["ytick.major.pad"] = 3.5

rcParams['legend.handlelength'] = 2.2

alpha = 1/137.036 # fine-structure constant

def m22_nu(m22):
    return 1e-6 * (m22/40)

def nu_m22(nu):
    return  40 * nu / 1e-6

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--arrows',action='store_true')
parser.add_argument('--summary',action='store_true')
parser.add_argument('--project',action='store_true')
args = parser.parse_args()


draw_arrows = args.arrows
draw_summary = args.summary
draw_project = args.project

filename = plot.get_datafile('fdm_limits.yaml')
limits = yaml.load(open(filename))

fig,ax = plt.subplots(figsize=(12,8))
ax.set_yscale('log'); ax.set_xscale('log')


# Cosmic Probes

if draw_summary:
    plot_limit_fill(limits['cosmic_probes'],zorder=10)
    plot_limit_fill(limits['cosmic_probes'],zorder=1,color='w',alpha=1)
    plot_limit(limits['cosmic_probes'],alpha=1.0,zorder=10)
    ax.annotate('Small-Scale Structure',(3e-21,1e-14),fontsize=16,rotation=90)
    ax.annotate('Extreme Environments',(1e-16,7e-12),fontsize=16)
    ax.annotate('CMB Polarization',(1e-22,5e-11),rotation=60,fontsize=16)
    ax.annotate('Stellar Interiors',(3e-2,6e-11),fontsize=16)
    ax.annotate('Decay Products',(5e1,1e-12),rotation=-75,fontsize=16)
    # BHSR
    plot_limit_fill(limits['bhsr_baryakhtar'],low=True,factor=alpha/(2*np.pi),alpha=0.5,zorder=1);
    plot_limit(limits['bhsr_baryakhtar'],factor=alpha/(2*np.pi),alpha=1.0,zorder=10);
    ax.annotate('Extreme Environments',(1e-12,7e-17),rotation=90,fontsize=16)
    #ax.annotate('Cosmic Probes',(1e-19,1e-12),color='g',fontsize=30)
    #ax.annotate('Cosmic Probes',(1e-18,1e-11),xycoords='data',fontsize=22)
    # CMB HD
    #plot_limit(limits['cmb_hd_oscillation'],color='r',alpha=1.0,zorder=10)
    #plot_limit(limits['cmb_hd_conversion'],color='r',alpha=1.0,zorder=10)
    #ax.annotate('CMB HD',(1e-17,7e-14),fontsize=16,color='r')
    ax.annotate(r'{\bf Cosmic Probes}',(1e-22,1e-9),color='k',fontsize=30,zorder=999)
else:
    #plot_limit_fill(limits['bhsr_mehta'],low=True,factor=alpha/(2*np.pi),alpha=0.5,zorder=1); plot_text(limits['bhsr_mehta'])
    plot_limit_fill(limits['bhsr_baryakhtar'],low=True,factor=alpha/(2*np.pi),alpha=0.5,zorder=1); plot_text(limits['bhsr_baryakhtar'])
    plot_limit_fill(limits['rogers2022_lya'],zorder=1); plot_text(limits['rogers2022_lya'])
    plot_limit_fill(limits['nadler2021_dsphs'],zorder=1); plot_text(limits['nadler2021_dsphs'])
    plot_limit_fill(limits['chen2022_m87'],zorder=1); plot_text(limits['chen2022_m87'])
    #plot_limit_fill(limits['xmm_newton'])
    plot_limit_fill(limits['ebl2'],zorder=1)
    plot_limit_fill(limits['leo_t'],zorder=1); plot_text(limits['leo_t'])
    plot_limit_fill(limits['xrays'],zorder=1); plot_text(limits['xrays'])
    plot_limit_fill(limits['ebl'],zorder=1); plot_text(limits['ebl'])
    plot_limit_fill(limits['chandra_ngc1275'],zorder=1); plot_text(limits['chandra_ngc1275'])
    plot_limit_fill(limits['mwd_polarization'],zorder=1); plot_text(limits['mwd_polarization'])
    #plot_limit_fill(limits['horizontal_branch'],zorder=1); plot_text(limits['horizontal_branch'])
    #plot_limit_fill(limits['horizontal_branch'],alpha=0.5,zorder=10)
    plot_limit_fill(limits['globular_clusters'],zorder=1); plot_text(limits['globular_clusters'])
    plot_limit_fill(limits['globular_clusters'],alpha=0.5,zorder=10)
    plot_limit_fill(limits['fedderke2019_cmb'],alpha=0.5,zorder=10); plot_text(limits['fedderke2019_cmb'])
    #plot_limit(limits['cosmic_probes'],alpha=1.0,color='r',zorder=999)
    #ax.annotate('Cosmic Probes',(1e-19,1e-14),color='g',fontsize=30)
    ax.annotate(r'{\bf Cosmic Probes}',(1e-22,1e-9),color='k',fontsize=30,zorder=999)
#ax.annotate(r'{\bf Cosmic Probes}',(1e-19,3e-14),color='g',fontsize=30)


if draw_arrows:
    arrow = FancyArrowPatch((2.1e-20, 1e-15), (2.5e-18, 1e-15), mutation_scale=70, fc='g', ec='none')
    ax.add_patch(arrow)
    #arrow = FancyArrowPatch((1e-11, 6e-12), (1e-11, 4e-13), mutation_scale=70, fc='g', ec='none')
    #ax.add_patch(arrow)
    arrow = FancyArrowPatch((1e-12, 7e-13), (1e-12, 4e-14), mutation_scale=70, fc='g', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((1.2e2, 3e-14), (1.2e0, 1e-14), mutation_scale=70, fc='g', ec='none')
    ax.add_patch(arrow)
elif draw_project:
    plot_limit(limits['cosmic_project'],zorder=10)
    plot_limit(limits['bhsr_project'],zorder=10)


# Haloscopes & Helioscopes
plot_limit_fill(limits['hscopes'],zorder=1)
plot_limit(limits['hscopes'],alpha=1.0,zorder=1)
ax.annotate('Helioscopes',(1e-10,3e-10),xycoords='data',fontsize=22)
ax.annotate('Haloscopes',(3e-6,1e-11),xycoords='data',fontsize=22,rotation=90)
if draw_arrows:
    arrow = FancyArrowPatch((1e-6, 1e-14), (1e-8, 1e-14), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((5e-5, 1e-12), (5e-3, 1e-12), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)
    arrow =FancyArrowPatch((1e-3, 7e-11), (1e-3, 7e-12), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)
    arrow =FancyArrowPatch((2e-5, 5e-15), (2e-5, 5e-16), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)


# QCD Axion
plot_limit(limits['ksvz'],zorder=0,lw=1.5)
plot_limit(limits['dfsz'],zorder=0,lw=1.5)
#plt.plot([1e-30,1e15],[1e-41,1e4],'--',color='0.5',lw=1.5,zorder=0)
#plt.plot([1e-30,1e15],[1e-39,1e6],'--',color='0.5',lw=1.5,zorder=0)
ax.annotate('QCD Axion',(1e-4,2e-13),xycoords='data',rotation=59,fontsize=18)

# ALP DM
plot_limit(limits['ringwald2013_alp_dm'],zorder=0)
mass,limit = get_mass_limit(limits['ringwald2013_alp_dm'])
interp = interp1d(np.log10(mass),np.log10(limit))
func = lambda m: 10**interp(np.log10(m))
x = np.array([1e-18,1e-11])
y = func(x)
yerr = y/1.5
line,cap,bar = ax.errorbar(x, y, yerr=yerr,uplims=True,lw=2.5,capsize=0,ls='none')
cap[0].set_markersize(10)


ax.annotate('ALP DM',(1e-16,1e-15),rotation=40,fontsize=18)

# Axes Limits
ax.set_xlim(1e-23,1e4)
ax.set_ylim(1e-20,1e-8)

ax.set_xticks([1e-22,1e-16,1e-10,1e-4,1e2])
ax.set_xticklabels([r'$10^{-22}$',r'$10^{-16}$',r'$10^{-10}$',r'$10^{-4}$',r'$10^{2}$'],fontsize=22)

ax.set_yticks([1e-18,1e-16,1e-14,1e-12,1e-10,1e-8])
ax.set_yticklabels([r'$10^{-18}$',r'$10^{-16}$',r'$10^{-14}$',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$'],fontsize=22)

ax.set_xlabel(r'Dark matter mass (eV/$c^2$)',fontsize=26,labelpad=8)
ax.set_ylabel(r'Axion-photon coupling, $g_{\phi\gamma}$ (GeV$^{-1}$)',fontsize=26,labelpad=12)

#ax.set_xlabel('$m_\phi$ (eV)',fontsize=26,labelpad=8)
#ax.set_ylabel('$g_{\phi\gamma}$ (GeV$^{-1}$)',fontsize=26,labelpad=12)

plt.subplots_adjust(top=0.90,bottom=0.12,left=0.15)

plt.tight_layout()

outfile = 'alp_photon_limits.pdf'
if draw_summary:
    plt.savefig(outfile.replace('.pdf','_summary.pdf'))
    plt.savefig(outfile.replace('.pdf','_summary.png'))
else:
    plt.savefig(outfile)
    plt.savefig(outfile.replace('.pdf','.png'))

#plt.show()
