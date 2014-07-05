"""
Created at 4:54 PM on 03 Jul 2014

Project: match
Subprogram: galsex.py

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''


'''

import os
import sys
import glob                                 # Find files

#Dashboard
#n = 10 # number of simulated galaxies to create per input image
sexdir = '/opt/local/bin/'                  # sextractor dir
galtab = 'galmodel.tab'
wdir = os.getcwd()+'/'


# Local Subdirectories
sexconfig = 'sexconfig/'
simcomb = 'simcomb/'
simgalset = 'simgalset/'
simcat = 'simcat/'
sexcat = 'sexcat/'


# Prefixes for output files
galcomb = 'galcomb'                                   # image with n sim gal added
galmodelset = 'galmodelset'
sim = 'sim'
sex = 'sex'         #


# Set number of itertions, sextractor config file, zeropoint magnitude, and verbosity
iteration = [1000, 10]
defaultsex = 'default.sex'
defaultsim = 'defaultsim.sex'
zp = 25.9463                                # zeropoint magnitude
output = sys.stdout
verbose = 1


# Scan directory for input fits files while excluding program files and alert if no .fits input images are present
images = glob.glob('*.fits')

l = len(images)


# Loop through science image versions
for w in xrange(l):

    # Loop through galcomb_w_x.fits images
    for x in xrange(iteration[0]):

        # Run sextractor in single image mode on new image (galcomb*.fits)
        cmnd = sexdir+'sex '+wdir+simcomb+galcomb+'_'+str(w)+'_'+str(x)+'.fits -c '+sexconfig+defaultsex
        cmnd += ' -mag_zeropoint '+str(zp)+' -catalog_name '
        cmnd += wdir+sexcat+sex+'_'+str(w)+'_'+str(x)+'.cat'
        os.system(cmnd)

        # Run sextractor in dual image mode:
        #   Find mag of galmodelset*.fits image (used for aperture correction)
        #   using the same apertures from the corresponding galcomb*.fits image
        cmnd2 = sexdir+'sex '+wdir+simcomb+galcomb+'_'+str(w)+'_'+str(x)+'.fits,'
        cmnd2 += wdir+simgalset+galmodelset+str(x)+'.fits -c '+sexconfig+defaultsim
        cmnd2 += ' -mag_zeropoint '+str(zp)+' -catalog_name '
        cmnd2 += wdir+simcat+sim+'_'+str(w)+'_'+str(x)+'.cat'
        os.system(cmnd2)

        if verbose:
            output.write('Generated sextractor catalog: '
                         + sex+'_'+str(w)+'_'+str(x)+'.cat \n')
