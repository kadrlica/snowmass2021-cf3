#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
__editor__ = "Ethan O. Nadler, Sidney Mau"
import os
import copy
import yaml
from collections import OrderedDict as odict

import numpy as np
from scipy.interpolate import interp1d
import pylab as plt

import matplotlib
import mpl_toolkits.axisartist as aa
from matplotlib.path import Path
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.patches import FancyArrowPatch

rcParams = matplotlib.rcParams
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
rcParams["xtick.direction"] = 'inout'
rcParams["ytick.direction"] = 'inout'
rcParams["xtick.major.size"] = 10
rcParams["ytick.major.size"] = 10
rcParams["xtick.minor.size"] = 5
rcParams["ytick.minor.size"] = 5
rcParams["xtick.major.width"] = 1.2
rcParams["ytick.major.width"] = 1.2
rcParams["xtick.minor.width"] = 1.0
rcParams["ytick.minor.width"] = 1.0
rcParams["xtick.major.pad"] = 7
rcParams["ytick.major.pad"] = 3.5
rcParams["axes.linewidth"] = 1.5
rcParams['legend.handlelength'] = 2.2


ticks = odict([
    [1e-21, 'zeV'],
    [1e-18, 'aeV'],
    [1e-15, 'feV'],
    [1e-12, 'peV'],
    [1e-9,  'neV'],
    [1e-6, r'$\mu$eV'],
    [1e-3,  'meV'],
    [1,     'eV'],
    [1e3,   'keV'],
    [1e6,   'MeV'],
    [1e9,   'GeV'],
    [1e12,  'TeV'],
    [1e15,  'PeV'],
    [1e20,  '$10 M_\odot$'],
])

fig,ax = plt.subplots(figsize=(14,2))
ax.set_xlim(-23,22)
ax.set_ylim(-0.1,1.1)
plt.subplots_adjust(left=0.01,right=0.99)

# Axes Limits

ax.set_xticks(np.log10(ticks.keys()))
ax.set_xticks(np.arange(-20,15),minor=True)
ax.set_xticklabels(ticks.values(),fontsize=22)
ax.set_xlabel(r'Dark Matter Mass',fontsize=26,labelpad=20)

# Turn off y-axis
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_yaxis().set_visible(False)

# Top labels
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_position(('data',0))
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()

a = .075
b = .075

kwargs = dict(color='k',zorder=10)
x0 = 17.5
x1 = 17.75
ax.plot([(x0+x1)/2.],[0], marker='s',color='w',markersize=6, zorder=5, clip_on=False)
ax.plot([x0-b,x0+b],[-a,+a], **kwargs)
ax.plot([x1-b,x1+b],[-a,+a], **kwargs)

#ax.plot(ax.get_xlim(),(0,0),'s',markersize=15,color='k')
plt.savefig('axis.png',bbox_inches='tight')
plt.savefig('axis.pdf',bbox_inches='tight')
