from __future__ import print_function
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
import cubehelix
'''
Do a tidal force calculation with the output from a MESA stellar model.
by Meredith Rawls, July 2015

**Input**
file from MESA with at least 4 columns as follows:
star_age (yrs), log_Teff (K), log_R (R_sun), cz_xm (mass of convective layer in M_sun)

2nd and 3rd file in slightly different format, for comparison:
star_age (yrs), M (total mass in M_sun), log_L (L_sun), log_Teff (K)

THIS IS DIFFERENT FROM THE DEFAULT MESA OUTFILE #sorrynotsorry
(YOU SHOULD DEFINITELY RUN tidalcalc_run1st.py FIRST)


**Output**
result of integral in the same way as reported in Verbunt & Phinney (1995)
plot for comparison with Verbunt & Phinney, Fig 1

'''
# ***** DEFINE AT MOST ONE OF THESE TWO THINGS ***** #
#Rnow = 8.33			# option to stop integration at some radius
idxstop = 340       # option to stop integration at some index
# ***** DEFINE AT MOST ONE OF THESE TWO THINGS ***** #

# Read in the MESA stellar model
infile = '../jeanshare/mesa_histshort_really.data'
times, logTeffs, logRs, cz_masses = np.loadtxt(infile, comments='#', dtype=np.float64, usecols=(0,1,2,3), unpack=True)
Teffs = np.power(10, logTeffs)
Rads = np.power(10, logRs)
#color_times = (np.log10(times) - np.log10(np.min(times))) / (np.log10(np.max(times)) - np.min(np.log10(times)))
cx = cubehelix.cmap(startHue=240,endHue=-300,minSat=1,maxSat=2.5,minLight=.3,maxLight=.8,gamma=.9)

# Read in other stellar models for comparison
compare1 = '../jeanshare/schaller_2.5Msun.txt'
compare2 = '../jeanshare/schaller_2.0Msun.txt'
const_sigma_sb = 5.6704e-5 # cgs units
convert_to_ergpersec = 3.839e33 # 1 solar luminosity in erg/sec
time1s, logT1s, logL1s, M1s = np.loadtxt(compare1, comments='#', dtype=np.float64, usecols=(0,3,2,1), unpack=True)
time2s, logT2s, logL2s, M2s = np.loadtxt(compare2, comments='#', dtype=np.float64, usecols=(0,3,2,1), unpack=True)
T1s = np.power(10, logT1s); T2s = np.power(10, logT2s)
L1s = np.power(10, logL1s); L2s = np.power(10, logL2s)
# calculate radii from Stefan-Boltzmann law, L = 4*pi*sigma*R^2*T^4
R1s = np.sqrt(L1s * convert_to_ergpersec / (4. * np.pi * const_sigma_sb * np.power(T1s,4.))) / 6.96e10
R2s = np.sqrt(L2s * convert_to_ergpersec / (4. *np.pi * const_sigma_sb * np.power(T2s,4.))) / 6.96e10
logR1s = np.log10(R1s)
logR2s = np.log10(R2s)

# Some hard-wired values for calculating delta-ln-e (not needed for just I)
Mtot = 2.16		# mass of other star in binary, assumed constant
f = 1.0			# fudge factors are fun (read the damn paper, it's ~1)
q = 0.99		# mass ratio
Porb = 171.277	# orbital period in days

def integrand(Teff, cz_mass, Rad):
	# see Equation 5 in Verbunt & Phinney 1995
	teff_part = np.power((Teff/4500.), (4./3.))
	mass_part = np.power(cz_mass, (2./3.))
	rad_part = np.power(Rad, 8.)
	return teff_part * mass_part * rad_part

def dlnecalc(f, Mtot, q, Int, Porb):
	# see Equation 6 in Verbunt & Phinney 1995
	return -1.7e-5 * f * np.power(Mtot, -11./3.) * q * np.power((1.+q), -5./3.) * Int * np.power(Porb, -16./3.)
	
# Evaluate the integral at each time point with Simpson's rule
# (this is different from just evaluating the inteGRAND at each point)
Intlist = []
Intlist1 = []
Intlist2 = []
for idx, time in enumerate(times):
	Intlist.append( simps(integrand(Teffs[0:idx+1], cz_masses[0:idx+1], Rads[0:idx+1]), x = times[0:idx+1]) )
for idx, time in enumerate(time1s):
	Intlist1.append( simps(integrand(T1s[0:idx+1], M1s[0:idx+1], R1s[0:idx+1]), x = time1s[0:idx+1]) )
for idx, time in enumerate(time2s):
	Intlist2.append( simps(integrand(T2s[0:idx+1], M2s[0:idx+1], R2s[0:idx+1]), x = time2s[0:idx+1]) )

# Tell python how many time steps to use in the integration
# (either by defining Rnow or idxstop earlier)
try:
    Rnow
    idxstop = np.argmax(Rads > Rnow) # saves the *first instance* that Rads exceeds Rnow
    print('Only integrating from time=0 up to time when R={0} (time={1:6e}, index={2})'.format(Rnow, times[idxstop], idxstop))
except:
    try:
        idxstop
        print('Only integrating from time=0 up to time={0:6e} (index={1}, R={2})'.format(times[idxstop], idxstop, Rads[idxstop]))
    except:
        idxstop = -1
        print('Integrating over all time')

# Final result of the integral
Int = Intlist[idxstop]

print('Value of I is {0}, so log I is {1}'.format(Int, np.log10(Int)))

# Evaluate the final expression
dlne = dlnecalc(f, Mtot, q, Int, Porb)

print('delta-ln-e is {0}'.format(dlne))
print('log[-(delta-ln-e)/f] = {0}'.format(np.log10(-1.*dlne/f)))

# First figure: I vs. logR
plt.figure(1)
plt.axis([-0.2,2.5,8,25])
plt.scatter(logRs[1:idxstop], np.log10(Intlist[1:idxstop]), marker='o', lw=0, c=times[1:idxstop], s=40, cmap=cx)#, label='MESA 2.16 Msun')
plt.plot(logR1s[1:], np.log10(Intlist1[1:]), marker='o', c='k', ls='none', mec='none', label='Schaller 2.5 Msun')
plt.plot(logR2s[1:], np.log10(Intlist2[1:]), marker='o', c='b', ls='none', mec='none', label='Schaller 2.0 Msun')
plt.xlabel('log R/Rsun')
plt.ylabel('log I')
plt.axvline(x=np.log10(8.33), c='k')
plt.axvline(x=np.log10(7.90), c='k')
plt.title('colorbar is time in years (times $10^8$)')
plt.legend(loc='upper left', numpoints=1, frameon=False)
plt.colorbar()

# Second figure: R vs. time
plt.figure(2)
plt.scatter(Rads[1:idxstop], times[1:idxstop], marker='o', lw=0, s=40, c=np.log10(Intlist[1:idxstop]), cmap=cx)
plt.axvline(x=8.33, c='k')
plt.axvline(x=7.90, c='k')
plt.xlabel('logR')
plt.ylabel('time')
plt.title('colorbar is the integral logI')
plt.colorbar()
plt.show()