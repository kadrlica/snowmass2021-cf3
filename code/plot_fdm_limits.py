#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import pylab as plt
import numpy as np
import yaml

from lsstplot import plot_limit, plot_limit_fill, plot_limit_patch, plot_lsst_limit

limits = yaml.load(open('data/fdm_limits.yaml'))

fig,ax = plt.subplots()
ax.set_yscale('linear')
ax.set_xscale('log')

plot_lsst_limit(limits['mw_sats'])

plot_limit(limits['nadler19'])
plot_limit(limits['armengaud17'])

plt.xlim(1e-22,1e-20)
plt.ylim(0.0,1.0)
plt.xlabel(r'Dark Matter Mass (eV)',fontsize=18)
plt.ylabel(r'Dark Matter Fraction',fontsize=18)
plt.subplots_adjust(top=0.95,bottom=0.12)

plt.savefig('fdm_limits.pdf')
plt.ion()

