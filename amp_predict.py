from __future__ import division
from __future__ import print_function
import numpy as np
'''
Use delta-nu, nu-max, and temp to estimate oscillation amplitudes.
Based on Corsaro et al. 2013.
'''
# Read in text file containing useful info
infile = '../../RGEB_seismicinfo.txt'

# Define columns
KICcol = 0
numaxcol = 1; numaxerrcol = 2
dnucol = 3; dnuerrcol = 4
teffcol = 5; tefferrcol = 6
ampcol = 7

# Read in data
KICs = np.loadtxt(infile, comments='#', usecols=(KICcol,), unpack=True)
numaxes, numaxerrs = np.loadtxt(infile, comments='#', usecols=(numaxcol,numaxerrcol), unpack=True)
dnus, dnuerrs = np.loadtxt(infile, comments='#', usecols=(dnucol,dnuerrcol), unpack=True)
teffs, tefferrs = np.loadtxt(infile, comments='#', usecols=(teffcol,tefferrcol), unpack=True)
amps = np.loadtxt(infile, comments='#', usecols=(ampcol,), unpack=True)

#for star in KICs:
#    print(star)

def ampguess(numax = 3100, dnu = 135, teff = 5777):
    '''
    Calculate amplitude given numax, dnu, and teff
    '''
    # Solar values for reference
    amp_s = 3.6
    numax_s = 3100
    dnu_s = 135
    teff_s = 5777
    # Constants from Corsaro et al. 2013 paper (M_{4,beta} specifically)
    s = 0.602; s_err = 0.008; r = 5.87; r_err = 0.14; t = 1.31; t_err = 0.02
    ampfrac = (numax/numax_s)**(2*s-3*t) * (dnu/dnu_s)**(4*t-4*s) * (teff/teff_s)**(5*s-1.5*t-r+0.2)
    amp = ampfrac*amp_s
    return amp
    

for star, numax, dnu, teff, amp in zip(KICs, numaxes, dnus, teffs, amps):
    ampcalc = ampguess(numax, dnu, teff)
    print('{0}: {1:.1f} predicted, {2:.1f} observed ({3:.0f} percent)'.format(int(star), ampcalc, amp, amp/ampcalc*100))
