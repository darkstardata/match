"""
Created on Sat Jan  4 19:42:20 2014

Project: match
Subprogram: match.py

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''
Matches detections in sextractor catalogs from
galcomb*.fits and galmodelset*.fits

Runtime ~ 3 hr
'''

import os
import sys
import numpy as np
from astropy.table import Table


#Dashboard
iteration = [5, 2000]
output = sys.stdout
verbose = 1
sexdir = '/opt/local/bin/'                  # sextractor dir

# Field and Filter selector
y = 0       # Field index
z = 1       # Filter index

fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
pixel_resolution = [0.06, 0.06, 0.03, 0.03, 0.03]
field = fields[y]
filt = filters[z]
pixres = pixel_resolution[z]

# Local Subdirectories and special input and output files
wdir = os.getcwd() + '/'
science = wdir + 'science/' + field + '/' + filt + '/'
simcomb = science + 'simcomb/'
simgal = science + 'simgal/'
simgalset = science + 'simgalset/'
simcat = science + 'simcat/'
sexcat = science + 'sexcat/'
sexmatch = science + 'sexmatch/'
sexconfig = science + 'sexconfig/'

outfile = 'sexmatch'
galmodel = 'galmodel'                                 # sim galaxy file prefix
sim = 'sim'
sex = 'sex'

# Read in sim gal data catalog
galtab = 'galmodel_' + field + '_' + filt + '.tab'
gal = Table.read(science + galtab, format='ascii.tab')

# Array for sim gal data
gid = gal[gal.colnames[0]]
gxpix = gal[gal.colnames[1]]
gypix = gal[gal.colnames[2]]
gmag = gal[gal.colnames[3]]
greff = gal[gal.colnames[4]]
gsersic = gal[gal.colnames[5]]
gba = gal[gal.colnames[6]]
gpa = gal[gal.colnames[7]]


# First Summing array for sextractor gal data
sxpix = []
sypix = []
smag = []
smagerr = []
sreff = []
sab = []
spa = []
sisoarea = []
sflags = []

ssimmag = []
ssimmagerr = []

# Loop through galcomb_w_x.fits images
for x in xrange(iteration[0]):

    # Read new sextractor catalog for galcomb images into table
    sc = Table.read(sexcat + sex + str(x) + '.cat', format='ascii.sextractor')

    xpix = sc[sc.colnames[1]]
    ypix = sc[sc.colnames[2]]
    mag = sc[sc.colnames[3]]
    magerr = sc[sc.colnames[4]]
    reff = sc[sc.colnames[5]]
    ab = sc[sc.colnames[6]]
    pa = sc[sc.colnames[7]]
    isoarea = sc[sc.colnames[8]]
    flags = sc[sc.colnames[9]]

    # Read in sextractor catalog for galmodelset images
    simc = Table.read(simcat + sim + str(x) + '.cat', format='ascii.sextractor')

    simmag = simc[simc.colnames[3]]
    #simmagerr = simc[simc.colnames[4]]

    # Define blank array in which to store sim gal and sex matches
    xmatch = np.zeros(iteration[1])
    ymatch = np.zeros(iteration[1])
    magmatch = np.zeros(iteration[1])
    magerrmatch = np.zeros(iteration[1])
    reffmatch = np.zeros(iteration[1])
    abmatch = np.zeros(iteration[1])
    pamatch = np.zeros(iteration[1])
    isoareamatch = np.zeros(iteration[1])
    flagsmatch = np.zeros(iteration[1])

    simmagmatch = np.zeros(iteration[1])
    #simmagerrmatch = np.zeros(iteration[1])

    # Match batch of sim gal with sextractor detections
    for z in xrange(iteration[1]):

        # Match x,y elements from simulated with sextractor detections
        xbin = np.asarray(np.where(np.logical_and(xpix >= gxpix[x * iteration[1] + z] - pixres,
                                                  xpix <= gxpix[x * iteration[1] + z] + pixres)))[0]

        ybin = np.asarray(np.where(np.logical_and(ypix >= gypix[x * iteration[1] + z] - pixres,
                                                  ypix <= gypix[x * iteration[1] + z] + pixres)))[0]


        # Set x and y positions found by sextractor to variables
        if len(np.intersect1d(xbin, ybin)) == 1:
            xmatch[z] = xpix[np.intersect1d(xbin, ybin)[0]]
            ymatch[z] = ypix[np.intersect1d(xbin, ybin)[0]]
            magmatch[z] = mag[np.intersect1d(xbin, ybin)[0]]
            magerrmatch[z] = magerr[np.intersect1d(xbin, ybin)[0]]
            reffmatch[z] = reff[np.intersect1d(xbin, ybin)[0]]
            abmatch[z] = ab[np.intersect1d(xbin, ybin)[0]]
            pamatch[z] = pa[np.intersect1d(xbin, ybin)[0]]
            isoareamatch[z] = isoarea[np.intersect1d(xbin, ybin)[0]]
            flagsmatch[z] = flags[np.intersect1d(xbin, ybin)[0]]

            simmagmatch[z] = simmag[np.intersect1d(xbin, ybin)[0]]
            #simmagerrmatch[z] = simmagerr[np.intersect1d(xbin, ybin)[0]]

        # Flag multiple detections
        if len(np.intersect1d(xbin, ybin)) > 1:
            xmatch[z] += -1
            ymatch[z] += -1
            magmatch[z] += -1
            magerrmatch[z] += -1
            reffmatch[z] += -1
            abmatch[z] += -1
            pamatch[z] += -1
            isoareamatch[z] += -1
            flagsmatch[z] += -1

            simmagmatch[z] += -1
            #simmagerrmatch[z] += -1

    # Append data to aggregate list
    sxpix.extend(xmatch)
    sypix.extend(ymatch)
    smag.extend(magmatch)
    smagerr.extend(magerrmatch)
    sreff.extend(reffmatch)
    sab.extend(abmatch)
    spa.extend(pamatch)
    sisoarea.extend(isoareamatch)
    sflags.extend(flagsmatch)

    ssimmag.extend(simmagmatch)
    #ssimmagerr.extend(simmagerrmatch)

    # Write data to incremental tables
    f = Table()
    f['x   '] = xmatch
    f['y   '] = ymatch
    f['simmag'] = simmagmatch
    f['mag '] = magmatch
    f['magerr'] = magerrmatch
    f['reff'] = reffmatch
    f['ab  '] = abmatch
    f['pa  '] = pamatch
    f['flags'] = flagsmatch

    # Format table data
    f['x   '].format = '%7.2f'
    f['y   '].format = '%7.2f'
    f['simmag'].format = '%5.2f'
    f['mag '].format = '%5.2f'
    f['magerr'].format = '%6.4f'
    f['reff'].format = '%4.2f'
    f['ab  '].format = '%7.2f'
    f['pa  '].format = '%5.2f'
    f['flags'].format = '%3d'

    f.write(sexmatch + outfile + str(x) + '.tab', format='ascii.tab')

    if verbose:
        output.write('\n data written to ' + outfile + str(x) + '.tab \n')


t = Table()
t['ID'] = gid

t['xpix'] = gxpix
t['ypix'] = gypix
t['mag'] = gmag
t['r_eff'] = greff
t['n_sersic'] = gsersic
t['b/a'] = gba
t['PA'] = gpa

t['sex_xpix'] = list(sxpix)
t['sex_ypix'] = list(sypix)
t['sex_mag'] = list(smag)
t['sex_mag-err'] = list(smagerr)

t['sim_mag'] = list(ssimmag)
#t['sim_mag-err'] = list(ssimmagerr)
t['sex_reff'] = list(sreff)
t['sex_a/b'] = list(sab)
t['sex_PA'] = list(spa)

t['isoarea'] = list(sisoarea)
t['sflags'] = list(sflags)

# Format table data
t['ID'].format = '%4d'

t['xpix'].format = '%7.2f'
t['ypix'].format = '%7.2f'
t['mag'].format = '%5.2f'
t['r_eff'].format = '%4.2f'
t['n_sersic'].format = '%5.2f'
t['b/a'].format = '%4.2f'
t['PA'].format = '%5.2f'

t['sex_xpix'].format = '%7.2f'
t['sex_ypix'].format = '%7.2f'
t['sex_mag'].format = '%5.2f'
t['sex_mag-err'].format = '%6.4f'
t['sim_mag'].format = '%5.2f'
#t['sim_mag-err'].format = '%6.4f'
t['sex_reff'].format = '%4.2f'
t['sex_a/b'].format = '%7.2f'
t['sex_PA'].format = '%5.2f'

t['isoarea'].format = '%3d'
t['sflags'].format = '%3d'



t.write(science + 'all' + outfile + '.tab', format='ascii.tab')

if verbose:
    output.write('\n data written to all' + outfile + '.tab \n')
