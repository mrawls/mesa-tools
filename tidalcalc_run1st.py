from __future__ import print_function
import numpy as np
import mesa as ms
'''
Reads a MESA history.data file
Prints out the columns you need for tidalcalc.py
NOTE: Only keeps timestamps where the star age is > 10^7 years

IMPORTANT: Manually define outfile, and MESA LOG directory below
(You can easily put the names of different columns if you want them instead)
'''

#outfile = '../../jeanshare/histshort_2.15v2.data'
#modeldir = ms.history_data('../../jeanshare/LOGS_2.15_v2')
outfile =   '../../RG_ELCmodeling/9291629/mesa/WORK_DIR/histshort_1.14.data'
modeldir =  '../../RG_ELCmodeling/9291629/mesa/WORK_DIR/LOGS_1.14'
modeldir = ms.history_data(modeldir)

ages = modeldir.get('star_age')
teffs = modeldir.get('log_Teff')
rads = modeldir.get('log_R')
czxms = modeldir.get('cz_xm')
#czbots = modeldir.get('cz_bot_mass')
#cztops = modeldir.get('cz_top_mass')
with open(outfile, 'w') as fout:
    for age, teff, rad, czxm in zip(ages, teffs, rads, czxms):
        if age > 1e7:
            print('{0:10e} {1:10e} {2:10e} {3:10e}'.format(age, teff, rad, czxm), file=fout)