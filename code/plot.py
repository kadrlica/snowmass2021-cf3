#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
__editor__ = "Ethan O. Nadler"

import os
from StringIO import StringIO
from collections import OrderedDict as odict
import numpy as np
import pylab as plt
import matplotlib

from matplotlib.patches import Ellipse, PathPatch
from matplotlib.path import Path
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'

COLORS = odict([
    ('blue','#1F77B4'),  # This is original T10 "tab:blue"
    ('gray','#7F7F7F'),  # This is original T10 "tab:gray"
    ('orange','orange'),
    ('gold','orange'),
    ('red','#FA0303'),   # T10 color is #D62728
    ('black','k'),   # T10 color is #D62728
])
ALPHA = 0.35
LINEWIDTH = 2

custom_blues = ["#99DDFF","#66CCFF", "#33BBFF", "#00AAFF", "#0088CC", "#006699", "#004466", "#000000"]

DEFAULTS = dict(color=COLORS['blue'],alpha=ALPHA,linewidth=LINEWIDTH)

def get_datadir():
    """Get the data directory. Could live here or up one level."""
    if os.path.exists('./data'):
        return './data'
    elif os.path.exists('../data'):
        return '../data'
    else:
        raise IOError("Data directory not found.")

def get_datafile(filename):
    """Get a data file."""
    datadir = get_datadir()
    path = os.path.join(datadir,filename)
    if not os.path.exists(path):
        raise IOError("Data file not found: %s"%path)
    return path

def setdefaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

def get_mass_limit(data, factor=1.0):
    mass,limit = np.genfromtxt(StringIO(data['xystring'])).T
    if data.get('mass_unit') == 'gram':
        mass /= 2e33
    elif data.get('mass_unit') == 'tev':
         mass *= 1e3
    limit *= factor
    return mass, limit

def plot_text(data,**kw):
    kwargs = dict(fontsize=data.get('fontsize',20), ha='center',  va='top',
                  rotation=data.get('rotation',0),zorder=100)
    kwargs.update(kw)
    return plt.text(float(data['label_x']),
                    float(data['label_y']),
                    data['style']['label'],
                    **kwargs)

def plot_limit(data, factor=1.0, **kw):
    mass,limit = get_mass_limit(data, factor=factor)
    kwargs = dict(**data['style'])
    kwargs.update(kw)
    kwargs.pop('label',None)
    plt.plot(mass, limit, **kwargs)

def plot_limit_fill(data, low=False, factor=1.0, **kw):
    mass,limit = get_mass_limit(data, factor=factor)

    kwargs = dict(**data['style'])
    kwargs.update(kw)
    setdefaults(kwargs,DEFAULTS)

    plt.fill_between(mass, limit, y2 = 1 if not low else 0,
                     edgecolor=kwargs.get('edgecolor'),
                     facecolor=kwargs['color'],
                     alpha=kwargs['alpha'],
                     zorder=kwargs.get('zorder'),
    )


def plot_limit_final(data, low=False, factor=1.0):
    kwargs = dict(**data['style'])
    setdefaults(kwargs,DEFAULTS)

    mass,limit = get_mass_limit(data, factor=factor)
    kwargs['alpha'] = 0.15
    plt.fill_between(mass, limit, y2 = 1 if not low else 0,
                     edgecolor=kwargs['color'],
                     facecolor=kwargs['color'],
                     alpha=kwargs['alpha'],
    )
    plot_text(data)

    kwargs['alpha'] = 1
    plt.plot(mass, limit, **kwargs)


def plot_limit_patch(data,**kw):
    kwargs = dict(**data['style'])
    kwargs.update(kw)
    setdefaults(kwargs,DEFAULTS)

    mass,limit = get_mass_limit(data)
    patch = PathPatch(Path(zip(mass, limit)),
                      edgecolor=kwargs.get('edgecolor','none'),
                      facecolor=kwargs['color'],
                      linestyle=kwargs['linestyle'],
                      alpha=kwargs['alpha'],
                      linewidth=kwargs['linewidth'],
                      zorder=kwargs.get('zorder')
                      )
    plt.gca().add_artist(patch)

def plot_one(data):
    kwargs = dict(**data['style'])
    setdefaults(kwargs,DEFAULTS)

    mass,limit = get_mass_limit(data)
    plt.fill_between(mass, limit, y2=1.0,
                     edgecolor=COLORS['blue'],
                     facecolor=COLORS['blue'],
                     alpha=kwargs['alpha'])
    plot_text(data)

def plot_two(data_loose, data_tight):
    kwargs = dict(**data_loose['style'])
    setdefaults(kwargs,DEFAULTS)

    mass_loose,limit_loose = get_mass_limit(data_loose)
    mass_tight,limit_tight = get_mass_limit(data_tight)

    # Interpolate the data to a common grid
    x_min = np.min((np.min(mass_loose), np.min(mass_tight)))
    x_max = np.max((np.max(mass_loose), np.max(mass_tight)))
    x = np.logspace(np.log10(x_min), np.log10(x_max), num=100)
    limit_loose_interp = np.interp(x, mass_loose, limit_loose)
    limit_tight_interp = np.interp(x, mass_tight, limit_tight)

    plt.fill_between(mass_loose, limit_loose, y2=1.0,
                     edgecolor=COLORS['blue'],
                     facecolor=COLORS['blue'],
                     alpha=kwargs['alpha'])
    plt.fill_between(x, limit_tight_interp, limit_loose_interp,
                     edgecolor='k',
                     linewidth=0,
                     facecolor='k',
                     alpha=0.07)
    plot_text(data_loose)
