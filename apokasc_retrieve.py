import numpy as np
from astropy.io import fits
'''
Quick scripts for interacting with APOKASC/APOGEE catalogs I downloaded
'''

KICs = [8702921., 9291629., 3955867., 10001167., 5786154., 7037405., 9246715., 9970396.]
TMIDs = ['2M19274322+3904194', '2M19315429+4232516', '2M19463960+4451110', 
         '2M19545035+4649589', '2M19074937+4656118']

#catfile = '../../APOKASC_cat_v3.1.2.txt'
#data = np.loadtxt(catfile, comments='#', usecols=(1,97,98,100,101,103,104,106,107,182,183))

#print(data[0])

#for line in data:
#    if line[0] in KICs:
#        print(line)

#fitsdata = '../../apogee-rc-DR12.fits'
fitsdata = '../../allStar+-dr12.fits'
hdulist = fits.open(fitsdata)
tbdata = hdulist[1].data
cols = hdulist[1].columns

#print(cols)

for idx, star in enumerate(tbdata['APOGEE_ID']):
    if star in TMIDs:
        print(star, tbdata['diso'][idx])

#fitsdata2 = '../../v601michael.fits'
#hdulist = fits.open(fitsdata2)
#tbdata = hdulist[1].data
#cols = hdulist[1].columns
#print(cols)

#for idx, star in enumerate(tbdata['APOGEE_ID']):
#    if star in TMIDs:
#        print(star, tbdata['DIST_SOL'][idx], tbdata['SIG_DISTSOL'][idx])
        
#print(tbdata['APOGEE_ID'], tbdata['DIST_SOL'], tbdata['SIG_DISTSOL'])

#DIST_SOL
#SIG_DISTSOL
#DIST_GCRADIAL
#SIG_GCRADIAL
#DIST_GCTOTAL
#SIG_GCTOTAL