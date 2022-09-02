"""
DEPRECATED: Plot showing the expected uncertainty for a positive
detection of PBHs with LSST/Roman
"""
import math
import yaml
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def get_sensitivity(fname="data/macho_limits.yaml"):
    """Load the MACHO limits from disc."""
    with open(fname, 'r') as ff:
        data = yaml.safe_load(ff)
        #Rubin
        rubin = np.array([float(xx) for xx in data["lsst_microlensing"]['xystring'].split()]).reshape(-1, 2)
        #Roman
        roman = np.array([float(xx) for xx in data["lsst_m31_microlensing"]['xystring'].split()]).reshape(-1, 2)
        return np.concatenate([roman, rubin])

def discovery(mass, fraction, LSSTfrac, label="$5-\sigma$ Discovery Constraints"):
    """Simple order of magnitude estimate.
    Map from number of events to detection fraction uses
    figure 2 of snowmass paper."""
    nexpectedevent = fraction/np.exp(LSSTfrac(np.log(mass)))
    #Poisson error bar for uncertainty
    frac = 5/np.sqrt(nexpectedevent)
    #X error is using the 50% error on the mass measurement of 2202.01903
    #plt.fill_between(mass, fraction *(1- frac), fraction *(1+frac), color="blue", alpha=0.5, edgecolor="white", label=label)
    #mm = 30
    #ff = 0.01
    yerror = fraction*np.array([frac, frac])
    xerror = np.array([mass - mass/1.5, mass*1.5-mass])
    #xerror = np.array([mm - mm/1.5, mm*1.5-mm])
    #print(yerror, frac2)
    subsample = 12
    plt.errorbar(mass[::subsample], fraction[::subsample], fmt='o', capsize=5, yerr=yerror[:,::subsample].reshape(2,-1), label=label, xerr=xerror[:,::subsample].reshape(2,-1))
    plt.plot(mass, fraction, label="PBH Population", color="black")
    plt.yscale('log')
    plt.xscale('log')

forecasts = get_sensitivity()
LSSTfrac = interp1d(np.log(forecasts[:,0]), np.log(forecasts[:,1]))
plt.fill_between(forecasts[:,0], forecasts[:, 1], np.ones_like(forecasts[:,0]), hatch="/", color="gold", facecolor="white", label="Forecast Rubin/Roman constraints")

masses = np.logspace(-3.0, 4, 100)

def pbh_mass_func(mass, mc = 30, sigma=0.5, fpbh = 0.05):
    """Lognormal mass function from 1705.05567 eq 3"""
    return fpbh / np.sqrt(2*np.pi) * np.exp(-1*np.log(mass/mc)**2/(2*sigma**2))


fractions = pbh_mass_func(masses)
discovery(masses, fractions, LSSTfrac)

#m2 = np.array([5e-11, 1e-2, 30])
#tfrac = truefrac * np.ones_like(m2)
#tfrac[0] = 0.02

#masses = np.logspace(-11, -9.6, 100)
#fractions = tfrac[0] * np.ones_like(masses)
#discovery(masses, fractions, LSSTfrac, label=None)

#plt.plot(m2, tfrac, 'o', markersize=12,color="black")

#discovery([30, 0.05)
#discovery(1e-9, 0.05, 1e-3)
#discovery(1e-3, 0.05, 1e-3)
plt.legend(loc="lower left")
plt.xlabel("Compact Object Mass ($M_\odot$)")
plt.ylabel("Dark Matter Fraction")
plt.ylim(1e-4, 1)
plt.xlim(1e-14, 1e4)
plt.savefig("pbh_discovery.pdf")

"""Discovery potential for three representative compact object masses. Shown is the expected uncertainty for a population of compact objects of a single mass making up $0.05$ of the dark matter. We assume Poisson uncertainties with the observational number counts expected from the relevant future microlensing experiment (Rubin or Roman). True compact object masses are $10^{-9} M_\odot$ (orange), $10^{-3} M_\odot$ (green), $30 M_\odot$ (blue). Uncertainty in the final mass is shown as $50\%$, conservatively using the uncertainty in the recent detection of a microlensing event by Ref.~\cite{2022arXiv220201903L}.

@ARTICLE{2022arXiv220201903L,
       author = {{Lam}, Casey Y. and {Lu}, Jessica R. and {Udalski}, Andrzej and {Bond}, Ian and {Bennett}, David P. and {Skowron}, Jan and {Mroz}, Przemek and {Poleski}, Radek and {Sumi}, Takahiro and {Szymanski}, Michal K. and {Kozlowski}, Szymon and {Pietrukowicz}, Pawel and {Soszynski}, Igor and {Ulaczyk}, Krzysztof and {Wyrzykowski}, Lukasz and {Miyazaki}, Shota and {Suzuki}, Daisuke and {Koshimoto}, Naoki and {Rattenbury}, Nicholas J. and {Hosek}, Matthew W., Jr. and {Abe}, Fumio and {Barry}, Richard and {Bhattacharya}, Aparna and {Fukui}, Akihiko and {Fujii}, Hirosane and {Hirao}, Yuki and {Itow}, Yoshitaka and {Kirikawa}, Rintaro and {Kondo}, Iona and {Matsubara}, Yutaka and {Matsumoto}, Sho and {Muraki}, Yasushi and {Olmschenk}, Greg and {Ranc}, Clement and {Okamura}, Arisa and {Satoh}, Yuki and {Ishitani Silva}, Stela and {Toda}, Taiga and {Tristram}, Paul J. and {Vandorou}, Aikaterini and {Yama}, Hibiki and {Abrams}, Natasha S. and {Agarwal}, Shrihan and {Rose}, Sam and {Terry}, Sean K.},
        title = "{An isolated mass gap black hole or neutron star detected with astrometric microlensing}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
         year = 2022,
        month = feb,
          eid = {arXiv:2202.01903},
        pages = {arXiv:2202.01903},
archivePrefix = {arXiv},
       eprint = {2202.01903},
 primaryClass = {astro-ph.GA},
       adsurl = {https://urldefense.proofpoint.com/v2/url?u=https-3A__ui.adsabs.harvard.edu_abs_2022arXiv220201903L&d=DwIGAg&c=gRgGjJ3BkIsb5y6s49QqsA&r=jqra9u-fMUhrjdZ82nfgRz_-t2bAIBXFoP0gVHU2oek&m=fslIJ7MxfFZh4lIH4SE4qjp_SCrCMUO5b6H09tx-HFLur8skLWQn8Zub3M4dn4iN&s=UCa542d4pdl4H7FU7C4-NG8C3tfBdSW9kpuSUKnoGxY&e= },
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
