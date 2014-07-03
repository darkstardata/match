"""
Created at 7:56 PM on 02 Jul 2014

Project: match
Subprogram: galsex

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''
Combines science images with galmodelset*.fits images and outputs galcomb_*_*.fitsfor processing with sextractor.

'''
import os
import sys
import glob                                 # Find files
import pyfits                               # Open fits files

#Dashboard
#n = 10 # number of simulated galaxies to create per input image
zp = 25.9463                                # zeropoint magnitude
pixres = 2                                  # ~2kpc x ~2kpc postage stamp
xx = 1089
yy = 963


# Local Subdirectories
simcomb = 'simcomb/'
simgalset = 'simgalset/'

# Prefixes for output files
galtab = 'galmodel.tab'
galcomb = 'galcomb'                                   # image with n sim gal added
galmodel = 'galmodel'                                 # sim galaxy file prefix
galmodelset = 'galmodelset'

outfile = 'sexmatch'
outfilename = 'sexmatch.tab'
wdir = os.getcwd()+'/'

# Set outfile name, number of itertions, and verbosity
iteration = [1000, 10]

# Scan directory for input fits files while excluding program files and alert if no .fits input images are present
images = glob.glob('*.fits')
if outfilename in images:
    images.remove(outfilename)

cut1 = glob.glob(galmodel+'*.fits')
rem1 = len(cut1)
for i in xrange(rem1):
    if cut1[i] in images:
        images.remove(cut1[i])

cut2 = glob.glob(galcomb+'*.fits')
rem2 = len(cut2)
for i in xrange(rem2):
    if cut2[i] in images:
        images.remove(cut2[i])

if len(images) == 0:
    print('No input images. Call with "--help" for help.\n')
    sys.exit()
l = len(images)


# Loop through science images
for w in xrange(l):

    # Loop through galmodelset images
    for x in xrange(iteration[0]):

        # Combine Science image with simulated batch image
        image = pyfits.open(images[w])[0].data
        gimages = pyfits.open(simgalset + galmodelset + str(x) + '.fits')[0].data
        sexinput = image+gimages

        # Make sure no file exists with the same name before writing
        cut3 = glob.glob(simcomb + galcomb + '*.fits')
        if simcomb + galcomb + str(x) + '.fits' in cut3:
            os.remove(simcomb + galcomb + '_' + str(w) + '_' + str(x) + '.fits')

        # Write science+simulated image to file
        pyfits.writeto(simcomb + galcomb + '_' + str(w) + '_' + str(x) + '.fits', sexinput)