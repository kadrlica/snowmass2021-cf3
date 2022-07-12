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

import pylab as plt
import matplotlib
from matplotlib.patches import FancyArrowPatch

import plot
from plot import plot_limit, plot_limit_fill, plot_limit_patch, get_mass_limit, plot_text
from plot import custom_blues

matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['legend.handlelength'] = 2.2

def m22_nu(m22):
    return 1e-6 * (m22/40)

def nu_m22(nu):
    return  40 * nu / 1e-6

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--arrows',action='store_true')
parser.add_argument('--summary',action='store_true')
args = parser.parse_args()


draw_arrows = args.arrows
draw_summary = args.summary

filename = plot.get_datafile('fdm_limits.yaml')
limits = yaml.load(open(filename))

fig,ax = plt.subplots(figsize=(12,8))
ax.set_yscale('log'); ax.set_xscale('log')


# Cosmic Probes

if draw_summary:
    plot_limit_fill(limits['cosmic_probes'],zorder=10)
    plot_limit(limits['cosmic_probes'],alpha=1.0,zorder=10)
    ax.annotate('Cosmic Probes',(1e-18,1e-11),xycoords='data',fontsize=22)
else:
    plot_limit_fill(limits['rogers2022_lya']); plot_text(limits['rogers2022_lya'])
    plot_limit_fill(limits['nadler2021_dsphs']); plot_text(limits['nadler2021_dsphs'])
    plot_limit_fill(limits['chen2022_m87']); plot_text(limits['chen2022_m87'])
    #plot_limit_fill(limits['xmm_newton'])
    plot_limit_fill(limits['ebl2'])
    plot_limit_fill(limits['leo_t']); plot_text(limits['leo_t'])
    plot_limit_fill(limits['xrays']); plot_text(limits['xrays'])
    plot_limit_fill(limits['ebl']); plot_text(limits['ebl'])
    plot_limit_fill(limits['mwd_polarization']); plot_text(limits['mwd_polarization'])
    plot_limit_fill(limits['horizontal_branch'],zorder=1); plot_text(limits['horizontal_branch'])
    plot_limit_fill(limits['horizontal_branch'],alpha=0.5,zorder=10)
    plot_limit_fill(limits['fedderke2019_cmb'],alpha=0.5,zorder=10); plot_text(limits['fedderke2019_cmb'])
    ax.annotate('Cosmic Probes',(1e-19,1e-12),color='g',fontsize=30)

if draw_arrows:
    arrow = FancyArrowPatch((2.1e-20, 1e-15), (2.5e-18, 1e-15), mutation_scale=50, fc='g', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((1e-11, 6e-12), (1e-11, 6e-13), mutation_scale=50, fc='g', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((1.2e2, 3e-14), (1.2e0, 1e-14), mutation_scale=50, fc='g', ec='none')
    ax.add_patch(arrow)


# Haloscopes & Helioscopes
plot_limit_fill(limits['hscopes'],zorder=1)
plot_limit(limits['hscopes'],alpha=1.0,zorder=1)
ax.annotate('Helioscopes',(1e-10,7e-10),xycoords='data',fontsize=22)
ax.annotate('Haloscopes',(3e-6,2e-11),xycoords='data',fontsize=22,rotation=90)
if draw_arrows:
    arrow = FancyArrowPatch((1e-6, 1e-14), (1e-8, 1e-14), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)
    arrow = FancyArrowPatch((5e-5, 1e-12), (5e-3, 1e-12), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)
    arrow =FancyArrowPatch((1e-3, 7e-11), (1e-3, 7e-12), mutation_scale=50, fc='0.7', ec='none')
    ax.add_patch(arrow)


# QCD Axion
plt.plot([1e-30,1e15],[1e-41,1e4],'k--',lw=1.5,zorder=0)
plt.plot([1e-30,1e15],[1e-39,1e6],'k--',lw=1.5,zorder=0)
ax.annotate('QCD Axion',(1e-4,3e-13),xycoords='data',rotation=59,fontsize=18)

# Axes Limits
ax.set_xlim(1e-23,1e4)
ax.set_ylim(1e-18,1e-8)

ax.set_xticks([1e-22,1e-16,1e-10,1e-4,1e2])
ax.set_xticklabels([r'$10^{-22}$',r'$10^{-16}$',r'$10^{-10}$',r'$10^{-4}$',r'$10^{2}$'],fontsize=22)

ax.set_yticks([1e-18,1e-16,1e-14,1e-12,1e-10,1e-8])
ax.set_yticklabels([r'$10^{-18}$',r'$10^{-16}$',r'$10^{-14}$',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$'],fontsize=22)

ax.set_xlabel('$m_\phi$ (eV)',fontsize=26,labelpad=8)
ax.set_ylabel('$g_{\phi\gamma}$ (GeV$^{-1}$)',fontsize=26,labelpad=12)

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
