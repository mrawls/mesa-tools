from __future__ import print_function
from __future__ import division
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
import cubehelix
import mesa as ms
'''
Do a tidal force calculation with the output from a MESA stellar model.
by Meredith Rawls, July 2015
update March 2016: you no longer need 'tidalcalc_run1st.py' because that was silly

**Input**
history outfile from MESA that contains the following 4 columns:
star_age (yrs), log_Teff (K), log_R (R_sun), cz_xm (mass of convective layer in M_sun)
NOTE cz_xm is NOT written out by default and you will need to turn it on!

**Output**
1) result of integral in the same way as reported in Verbunt & Phinney (1995)
2) plot for comparison with Verbunt & Phinney, Fig 1

'''
# SET THESE VALUES CORRECTLY!
#sysname = '9246715'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_2.15_cz'
#modeldir = '../../jeanshare/LOGS_2.15_v2'
#Mstar = 2.15        # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 8.30        # present radius for plot
#Mcomp = 2.17        # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass =  0.01385255858
#Porb = 171.27688    # orbital period in days, needed for delta-ln-e
#t0 = 1e7            # age in years to begin integrating star (set near end of MS phase)

sysname = '8702921'
modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.59'
Mstar = 1.59        # mass of main star used in MESA model, needed for delta-ln-e
Rstar = 5.24        # present radius for plot
Mcomp = 0.272       # mass of other star in binary, assumed constant, needed for delta-ln-e
Zmass = 0.017100
Porb = 19.38446     # orbital period in days, needed for delta-ln-e
t0 = 1.7e9          # age in years to begin integrating star (set near end of MS phase)

#sysname = '9291629'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.1274'
#Mstar = 1.1274      # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 7.932       # present radius for plot
#Mcomp = 1.1147      # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass = 0.014658
#Porb = 20.686424    # orbital period in days, needed for delta-ln-e
#t0 = 5e9            # age in years to begin integrating star (set near end of MS phase)

#sysname = '3955867'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.103'
#Mstar = 1.103       # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 8.238       # present radius for plot
#Mcomp = 0.9187      # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass = 0.003475
#Porb = 33.65685     # orbital period in days, needed for delta-ln-e
#t0 = 7e7            # age in years to begin integrating star (set near end of MS phase)

#sysname = '10001167'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_0.857'
#Mstar = 0.857       # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 12.73       # present radius for plot
#Mcomp = 0.805       # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass = 0.002532
#Porb = 120.3903     # orbital period in days, needed for delta-ln-e
#t0 = 2e9          # age in years to begin integrating star (set near end of MS phase)

#sysname = '5786154'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.062'
#Mstar = 1.062       # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 11.01       # present radius for plot
#Mcomp = 1.019       # mass of other star in binary, assumed constant, needed for delta-ln-e
#Porb = 197.9182     # orbital period in days, needed for delta-ln-e
#Zmass = 0.010637
#t0 = 5e7          # age in years to begin integrating star (set near end of MS phase)

#sysname = '7037405'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.267'
#Mstar = 1.267       # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 13.72       # present radius for plot
#Mcomp = 1.15        # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass = 0.005689
#Porb = 207.1082     # orbital period in days, needed for delta-ln-e
#t0 = 2e9          # age in years to begin integrating star (set near end of MS phase)

#sysname = '9970396'
#modeldir =  '../../RG_ELCmodeling/'+sysname+'/mesa/LOGS_1.17'
#Mstar = 1.17        # mass of main star used in MESA model, needed for delta-ln-e
#Rstar = 7.70        # present radius for plot
#Mcomp = 0.998       # mass of other star in binary, assumed constant, needed for delta-ln-e
#Zmass = 0.007167
#Porb = 235.29852    # orbital period in days, needed for delta-ln-e
#t0 = 2e9          # age in years to begin integrating star (set near end of MS phase)

f = 1.0	        # fudge factors are fun (read the paper, it's ~1), needed for delta-ln-e
idxstart = 10   # timestamp at which to begin plotting points (default = 20)
q = Mcomp/Mstar # mass ratio = M_comp/M

# ***** DEFINE AT MOST ONE OF THESE TWO THINGS ***** #
Rnow = Rstar  		# option to stop integration at some radius
#idxstop = 2890      # option to stop integration at some index
# ***** DEFINE AT MOST ONE OF THESE TWO THINGS ***** #

# Read in MESA stellar model data from a history file
# Save only values corresponding to age > t0 years
model = ms.history_data(modeldir)
times = model.get('star_age')
t0idx = np.argmax(times > t0)
times = times[t0idx:]
logTeffs = model.get('log_Teff')[t0idx:]
logRs = model.get('log_R')[t0idx:]
cz_masses = model.get('cz_xm')[t0idx:]
loggs = model.get('log_g')[t0idx:]

Teffs = np.power(10, logTeffs)
Rads = np.power(10, logRs)


#len1 = len(times) #original length of timestamps
#times = np.array([time for time in times if time >= 1e7])
#len2 = len(times) #truncated length of timestamps >= 1e7
#Teffs = np.array([np.power(10,logTeff) for idx, logTeff in enumerate(logTeffs) if idx >= len1-len2])
#Rads = np.array([np.power(10,logR) for idx, logR in enumerate(logRs) if idx >= len1-len2])
#cz_masses = np.array([cz_mass for idx, cz_mass in enumerate(cz_masses) if idx >= len1-len2])

# colorbar settings
cx = cubehelix.cmap(startHue=240,endHue=-300,minSat=1,maxSat=2.5,minLight=.3,maxLight=.8,gamma=.9)

def integrand(Teff, cz_mass, Rad):
	'''
	see Equation 5 in Verbunt & Phinney 1995
	'''
	teff_part = np.power((Teff/4500.), (4./3.))
	mass_part = np.power(cz_mass, (2./3.))
	rad_part = np.power(Rad, 8.)
	return teff_part * mass_part * rad_part

def dlnecalc(f, Mcomp, q, Int, Porb):
	'''
	see Equation 6 in Verbunt & Phinney 1995
	'''
	return -1.7e-5 * f * np.power(Mcomp, -11./3.) * q * np.power((1.+q), -5./3.) * Int * np.power(Porb, -16./3.)
	
# Evaluate the integral at each time point with Simpson's rule
# (this is different from just evaluating the inteGRAND at each point, d'oh)
Intlist = []
for idx, time in enumerate(times):
	Intlist.append( simps(integrand(Teffs[0:idx+1], cz_masses[0:idx+1], Rads[0:idx+1]), x = times[0:idx+1]) )

# Tell python how many time steps to use in the integration
# (either by defining Rnow or idxstop earlier)
try:
    Rnow    
except:
    try:
        idxstop
    except:
        idxstop = -2
        print('Integrating over all time')
    else:
        print('Only integrating from time={0:6e} up to time={1:6e} (index={2}, R={3})'.format(times[0], times[idxstop], idxstop, Rads[idxstop]))
else:
    idxstop = np.argmax(Rads > Rnow) # saves the *first instance* that Rads exceeds Rnow
    print('Only integrating from time={0:6e} up to time when R={1} (time={2:6e}, index={3})'.format(times[0], Rnow, times[idxstop], idxstop))

# Final result of the integral
Int = Intlist[idxstop]
print('Value of I is {0}, so log I is {1}'.format(Int, np.log10(Int)))
print('Teff = {0}, logg = {1}'.format(Teffs[idxstop], loggs[idxstop]))

# Evaluate the final expression
dlne = dlnecalc(f, Mcomp, q, Int, Porb)

print('delta-ln-e is {0}'.format(dlne))
print('log[-(delta-ln-e)/f] = {0}'.format(np.log10(-1.*dlne/f)))

# Plot setup
fig = plt.figure()
fig.text(0.100, 0.93, 'KIC {0}'.format(sysname), size=30)
fig.text(0.373, 0.93, '$P_{{orb}}$ = {0:.1f} d'.format(Porb), size=24)
fig.text(0.573, 0.93, '$M$ = {0:.2f} $M_{{\odot}}$'.format(Mstar), size=24)
fig.text(0.773, 0.93, '$Z$ = {0:.4f}'.format(Zmass), size=24)
#fig.text(0.773, 0.93, '$M_2/M$ = {0:.3f}'.format(q), size=24)

# Top pane: R vs. time
ax2 = fig.add_subplot(2,1,1)
plt.scatter(logRs[idxstart:idxstop], times[idxstart:idxstop], 
            marker='o', lw=0, s=40, c=np.log10(Intlist[idxstart:idxstop]), cmap=cx)
plt.axvline(x=np.log10(Rstar), c='k', label='$R$ = {0} $R_{{\odot}}$'.format(Rstar))
plt.ylabel('Age (yr)', size=30)
color2 = plt.colorbar()
color2.set_label('log I (yr)', size=30)
#plt.text(1.2, 0.0, 'vertical line R = {0} Rsun'.format(Rstar))
plt.legend(loc=2, frameon=False) #loc=5 is center right, 10 is center, 2 is upper left
ax2.set_xticklabels([])

# Bottom pane: I vs. logR
ax1 = fig.add_subplot(2,1,2)
plt.subplots_adjust(hspace=0.1)
plt.scatter(logRs[idxstart:idxstop], np.log10(Intlist[idxstart:idxstop]), 
            marker='o', lw=0, c=times[idxstart:idxstop], s=40, cmap=cx)
plt.ylabel('log I (yr)', size=30)
plt.axvline(x=np.log10(Rstar), c='k')#, label='$R$ = {0} $R_{{\odot}}$'.format(Rstar))
plt.xlabel('log $R/R_{{\odot}}$', size=30)
color1 = plt.colorbar()
color1.set_label('Age (yr)', size=30)

plt.show()