#!/usr/bin/env python
"""
Plot dark matter mass vs indirect detection annihilation cross section.
"""
import pylab as plt
import numpy as np
import yaml

from lsstplot import plot_one, plot_two, plot_limit, plot_lsst_limit
from plot import plot_limit, plot_limit_fill, plot_limit_patch, get_mass_limit, plot_text


import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--arrows',action='store_true')
parser.add_argument('--summary',action='store_true')
parser.add_argument('--project',action='store_true')
args = parser.parse_args()

fig,ax = plt.subplots()
ax.set_yscale('log')
ax.set_xscale('log')

limits = yaml.load(open('data/macho_limits.yaml'))


# plot_lsst_limit(limits['lsst_microlensing'])
# plot_limit(limits['lsst_paralensing'])

if args.summary:
    # Aggressive (tight) limits
    plot_limit_fill(limits['cosmic_macho_tight'],color='0.8',zorder=0)
    plot_limit_fill(limits['hsc_niikura_2017'],color='0.8',zorder=0)
    # Conservative (loose) limits
    plot_limit_fill(limits['gammaray_background_loose_carr_2016'],color='g')
    plot_limit(limits['gammaray_background_loose_carr_2016'],color='g',linewidth=3)

    plot_limit_fill(limits['cosmic_macho_loose'],color='w',zorder=1,alpha=1)
    plot_limit_fill(limits['cosmic_macho_loose'],color='g',zorder=2)
    plot_limit(limits['cosmic_macho_loose'],color='g',linewidth=3,alpha=1.0)

    if True:
        # Full text
        ax.annotate("Gamma-ray Background",(2e-18,3e-1),color='k',fontsize=14,va='top',rotation=86)
        ax.annotate("Milky Way"+"\n"+"Microlensing",(7e-3,2e-1),color='k',fontsize=14,va='top',ha='center')
        ax.annotate("Andromeda"+"\n"+"Microlensing",(7e-9,2e-1),color='k',fontsize=14,va='top',ha='center')

        ax.annotate("CMB",(1e4,3e-4),color='k',fontsize=14,rotation=-85)
        ax.annotate("Dwarf"+"\n"+"Galaxies",(2e4,5e-3),color='k',fontsize=14,ha='center',rotation=-85)
        ax.annotate("Wide Binaries",(3e5,5e-1),color='k',fontsize=14,ha='center')
        ax.annotate("Supernova"+"\n"+"Lensing",(1e4,2e-1),color='k',fontsize=14,ha='center')
        ax.annotate("Disk"+"\n"+"Stability",(5e6,1e-2),color='k',fontsize=14,ha='center')
        ax.annotate("X-Ray"+"\n"+"Binaries",(3e6,1e-3),color='k',fontsize=14,ha='center')
    else:
        # Abbreviations
        ax.annotate("GB",(2e-18,1e-2),color='k',fontsize=14,va='top',rotation=86)
        ax.annotate("MWML",(1e-3,1e-1),color='k',fontsize=14,va='top')
        ax.annotate("M31ML",(3e-10,3e-2),color='k',fontsize=14,va='top')

        kw = dict(ha='center',  va='top')
        ax.annotate("CMB",(1e4,2e-4),color='k',fontsize=14,rotation=-85)
        ax.annotate("DH",(1e4,6e-4),color='k',fontsize=14,**kw)
        ax.annotate("XB",(1e5,5e-4),color='k',fontsize=14,**kw)
        ax.annotate("EII",(5e3,4e-2),color='k',fontsize=14,**kw)
        ax.annotate("DS",(4e7,2e-2),color='k',fontsize=14,**kw)
        ax.annotate("WB",(1e5,6e-1),color='k',fontsize=14,**kw)
        ax.annotate("LSN",(1e3,5e-1),color='k',fontsize=14,**kw)


else:
    plot_lsst_limit(limits['lsst_m31_microlensing'])
    # plot_limit(limits['lsst_paralensing'])
    plot_lsst_limit(limits['lsst_microlensing'])
    plot_lsst_limit(limits['gammaray_forecast_coogan2020'])

    plot_two(limits['gammaray_background_loose_carr_2016'],
             limits['gammaray_background_tight_carr_2016'])
    # plot_one(limits['gammaray_femtolens_carr_2016'])
    # plot_two(limits['ns_capture_loose_capela_2013'],
    #          limits['ns_capture_tight_capela_2013'])
    plot_two(limits['hsc_smyth_2020'],
             limits['hsc_niikura_2017'])
    plot_two(limits['eridanus_brandt_2016_loose'],
             limits['eridanus_brandt_2016_tight'])
    plot_two(limits['lensed_sn_garcia-bellido_2017_loose'],
             limits['lensed_sn_zumalacarregui_2018_tight'])
    plot_two(limits['lensed_sn_garcia-bellido_2017_loose'],
             limits['lensed_sn_zumalacarregui_2018_tight'])
    plot_one(limits['eros_macho_blaineau_2022'])
    plot_two(limits['xraybinary_2017_loose'],
             limits['xraybinary_2017_tight'])
    plot_two(limits['dwarfheating_2021_loose'],
             limits['dwarfheating_2021_tight'])
    plot_two(limits['binaries_quinn_2009_loose'],
             limits['binaries_yoo_2003_tight'])
    plot_two(limits['plank_ali-haimoud_2016_loose'],
             limits['cmb_ricotti_2008_tight'])
    plot_one(limits['disk_lacey_1985'])

if args.project:
    plot_limit(limits['gammaray_forecast_coogan2020'],linestyle='--')
    #plot_text(limits['gammaray_forecast_coogan2020'],fontsize=16)
    ax.annotate("MeV Telescopes",(1e-16,2e-3),color='k',fontsize=14,va='top',rotation=86)

    plot_limit(limits['lsst_m31_microlensing'],linestyle='--')
    ax.annotate("Roman Microlensing",(1e-12,5e-4),color='k',fontsize=14,va='top')

    plot_limit(limits['lsst_microlensing'],linestyle='--')
    ax.annotate("Rubin Microlensing",(1e-2,3.2e-4),color='k',fontsize=14,va='top')


plt.xlim(1e-18, 1e8)
plt.ylim(1e-4, 1.0)
plt.xlabel(r'${\rm Compact\ Object\ Mass}\ (M_\odot)$',fontsize=18)
plt.ylabel(r'${\rm Dark\ Matter\ Fraction}$', fontsize=18)
plt.subplots_adjust(top=0.95, bottom=0.12)

outbase = 'macho_limits'
if args.summary:
    outfile = outbase + '_summary.pdf'
else:
    outfile = outbase + '.pdf'

plt.savefig(outfile,rasterize=True)
plt.savefig(outfile.replace('.pdf','.png'),rasterize=True)


plt.ion()
