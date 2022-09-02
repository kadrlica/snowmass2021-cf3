"""
Plot the expected signal for a positive detection of PBHs with LSST.

References:
1709.07467 - Bellomo "Primordial Black Holes as Dark Matter: Converting Constraints from Monochromatic to Extended Mass Distributions"
1705.05567 - Carr "Primordial black hole constraints for extended mass functions"

Discussion on GitHub:
https://github.com/lsstdarkmatter/dark-matter-paper/issues/8
"""

import math
import yaml
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats
from scipy.integrate import quad
from scipy.interpolate import interp1d

from lsstplot import plot_one, plot_two, plot_limit, plot_lsst_limit

def plot_macho_limits(filename="macho_limits.yaml"):
    """ Plot existing and projected limits.

    Parameters
    ----------
    fiilename : filename

    Returns
    -------
    None
    """
    limits = yaml.load(open(filename))

    plot_lsst_limit(limits['lsst_microlensing'])

    plot_two(limits['gammaray_background_loose_carr_2016'],
             limits['gammaray_background_tight_carr_2016'])
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
    #plot_one(limits['disk_lacey_1985'])



def number_of_events(f_pbh):
    """Find expected number of events from PBHs.
    NOTE this neglects any potential dependence on optical depth
    from PBH mass for large PBHs which have an einstein radius
    large enough so that the crossing time is longer than survey time.

    https://github.com/lsstdarkmatter/dark-matter-paper/issues/8
    """
    #This is a lower bound (Will Dawson, private comm. :) )
    number_of_stars = 1e9
    #Polchinski 96 says 5e-7 to the LMC
    #We use Sumi 2005 et al (Ogle)
    optical_depth = 4.48e-6
    total = optical_depth * number_of_stars * f_pbh
    return total

def get_sensitivity(filename="macho_limits.yaml"):
    """Load the MACHO limits.

    Parameters
    ----------
    filename : yaml file with sensitivity

    Returns
    -------
    arr : array of mass,sensitivity
    """
    with open(filename, 'r') as f:
        data = yaml.safe_load(f)
        #Rubin
        rubin = np.array([float(x) for x in data["lsst_microlensing"]['xystring'].split()]).reshape(-1, 2)
        #Roman
        roman = np.array([float(x) for x in data["lsst_m31_microlensing"]['xystring'].split()]).reshape(-1, 2)
        return np.concatenate([roman, rubin])

def pbh_mass_func(mass, mc = 30., sigma=0.5, fpbh = 0.05):
    """Differential log-normal mass function from Eq 3 of 1705.05567.
    When integrated over all masses, this should give fpbh.

    Parameters
    ----------
    mass : PBH mass to evaluate the function
    mc   : characteristic mass
    sigma: characteristic width
    fpbh : normalization

    Returns
    -------
    dfpbh_dm : differential mass function
    """
    return fpbh * (1/np.sqrt(2*np.pi)) * 1/(sigma * mass) * np.exp(-1*np.log(mass/mc)**2/(2*sigma**2) )

def integrate_pbh_mass_func(mass_bins, mass_func):
    """

    Parameters
    ----------
    mass_bins : Bins of PBH mass (Msun)
    mass_func : Differential mass functions

    Returns
    -------
    mass : Centers of mass bins
    fpbh : Integrated PBH fraction
    """
    mass_centers = (mass_bins[1:] + mass_bins[:-1])/2.

    fpbh = []
    for i, (mlo, mhi) in enumerate(zip(mass_bins[:-1],mass_bins[1:])):
        fpbh.append(quad(mass_func,mlo,mhi)[0])
    return mass_centers, np.array(fpbh)


def pbh_expected_events(mass_bins, mass_func, lsst_sensitivity):
    """
    Calculate expected fpbh and number of PBH events.

    Parameters
    ----------
    mass_bins   : Edges of PBH mass bins
    mass_func   : PBH mass function, d(fpbh)/dM
    sensitivity : LSST sensitivity as a fraction of DM in PBHs

    Returns
    -------
    masses : central masses
    fpbh   : integrated fpbh
    nobs   : expected number of observed PBHs
    """

    masses, fpbh = integrate_pbh_mass_func(mass_bins, mass_func)

    # LSST sensitivity corresponds to measuring 1 PBH at that mass
    fpbh_limit = lsst_sensitivity(masses)

    # Number of expected PBHs expected at each bin
    nobs = fpbh / fpbh_limit

    return masses, fpbh, nobs

filename='data/macho_limits.yaml'
forecasts = get_sensitivity(filename)

def lsst_sensitivity(mass):
    """ Interpolation of LSST sensitivity

    Parameters
    ----------
    mass : PBH mass

    Returns
    -------
    fpbh : LSST fpbh sensitivity
    """
    loginterp = interp1d(np.log(forecasts[:,0]), np.log(forecasts[:,1]))
    return np.exp( loginterp(np.log(mass)) )


# Set up the mass function parameters
MC = 30 # Msun
SIGMA = 1.0 # width
FPBH = 0.03 # fraction
mass_func = lambda mass: pbh_mass_func(mass, MC, SIGMA, FPBH)

# Check that the pbh mass function is properly normalized
assert np.allclose(quad(mass_func,0,np.inf)[0], FPBH)

print("PBH mass function parameters:")
print("  fpbh:  %.3f"%FPBH)
print("  Mc   : %.1f Msun"%MC)
print("  sigma: %.1f"%SIGMA)

ntotal = number_of_events(FPBH)
print("PBH events: %.0f"%ntotal)

## Integrate fpbh more finely
full_masses = np.logspace(-1,4,100)
#full_masses, full_fpbh = integrate_pbh_mass_func(np.logspace(-1,4,100), mass_func)

# Calculate number of expected PBH events and fpbh in each mass bin
mass_bins = np.logspace(0, 3, 11)
masses, fpbh, nobs = pbh_expected_events(mass_bins, mass_func, lsst_sensitivity)
nobs_err = np.sqrt(nobs)
fpbh_err = fpbh * np.sqrt(nobs)/nobs
mass_err = [masses - mass_bins[:-1], mass_bins[1:] - masses]

uplims = nobs <= 1
nobs_err[uplims] = 1.0
fpbh_err[uplims] = 1.0

print("PBH events observed by LSST: %.0f"%(nobs.sum()))

color='k'

fig,ax = plt.subplots(1,1,figsize=(8,5))
plt.subplots_adjust(top=0.95,bottom=0.15)
plt.errorbar(masses,nobs,yerr=nobs_err,xerr=mass_err,uplims=uplims,capsize=5);
plt.axvline(MC,ls='--',color='gray')
plt.ylim(-2,None)
plt.xscale('log')
plt.xlabel(r"Compact Object Mass ($M_\odot$)")
plt.ylabel(r"Number of Compact Objects")
plt.savefig("pbh_nobserved.pdf")

fig,ax = plt.subplots(1,1,figsize=(8,5))
ax.set_yscale('log')
ax.set_xscale('log')

plot_macho_limits("data/macho_limits.yaml")

plt.errorbar(masses[~uplims], fpbh[~uplims], yerr=fpbh_err[~uplims],
             xerr=[mass_err[0][~uplims],mass_err[1][~uplims]],
             fmt='o', capsize=5,  color=color)

plt.xlim(1e-6, 1e6)
plt.ylim(1e-4, 1.0)
plt.xlabel(r'${\rm Compact\ Object\ Mass}\ (M_\odot)$',fontsize=18)
plt.ylabel(r'${\rm Dark\ Matter\ Fraction}$', fontsize=18)
plt.subplots_adjust(top=0.95, bottom=0.12)

plt.savefig("pbh_discovery.pdf")

