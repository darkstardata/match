# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 19:42:20 2014

Author: Andrew Crooks ; Graduate Student @ UC Riverside


"""

import os
import sys
import glob                                 # Find files
import numpy as np
import pyfits                               # Open fits files
from optparse import OptionParser           # Parse options
from astropy.table import Table
from astropy.io import fits

#Dashboard
#n = 10 # number of simulated galaxies to create per input image
zp = 25.9463                                # zeropoint magnitude
pixres = 2                                  # ~2kpc x ~2kpc postage stamp
xx = 1089
yy = 963

sexdir = '/opt/local/bin/'                  # sextractor dir
galmodel = 'galmodel'                       # sim galaxy file prefix
galtab = 'galmodel.tab'

# Local Subdirectories
simimages = 'simimages/'
simgalsets = 'simgalsets/'
simcat = 'simcat/'
sexcat = 'sexcat/'
sexmatch = 'sexmatch/'

# Prefixes for output files
mc = 'mc'                                   # image with n sim gal added
gsum = 'gsum'
sim = 'sim'
sex = 'sex'         #
outfile = 'sexmatch'
wdir = os.getcwd()+'/'


# Set outfile name, number of iterations, and verbosity
outfilename = 'sexmatch.tab'
iter = [1000, 10]
sexfile = 'default.sex''
verbose = 1
output = sys.stdout


# Scan directory for input fits files while excluding program files
images = glob.glob('*.fits')
if outfilename in images:
    images.remove(outfilename)

cut1 = glob.glob(galmodel+'*.fits')
rem1 = len(cut1)
for i in xrange(rem1):
    if cut1[i] in images:
        images.remove(cut1[i])

cut2 = glob.glob(mc+'*.fits')
rem2 = len(cut2)
for i in xrange(rem2):
    if cut2[i] in images:
        images.remove(cut2[i])

if len(images) == 0:
    print('No input images. Call with "--help" for help.\n')
    sys.exit()
l = len(images)


def main():

    # Read in sim gal data catalog
    gal = Table.read(galtab, format='ascii.tab')
    
    # Array for sim gal data
    gid = gal[gal.colnames[0]]
    gxpix = gal[gal.colnames[1]]
    gypix = gal[gal.colnames[2]]
    gmag = gal[gal.colnames[3]]
    greff = gal[gal.colnames[4]]
    gellipt = gal[gal.colnames[5]]
    gsersic = gal[gal.colnames[6]]
    gab = np.ones(len(gellipt)) - gellipt
    gpa = gal[gal.colnames[6]]

    # Loop through science images
    for w in xrange(l):

        # First Summing array for sextractor gal data
        sxpix = []
        sypix = []
        smag = []
        smagerr = []
        sflux = []
        sellipt = []
        sflags = []
        sisoarea = []

        ssimmag = []
        #ssimmagerr = []
        #ssimflux = []

        # Loop through galmodeset images
        for x in xrange(iter[0]):
            image = pyfits.open(images[w])[0].data
            gimages = pyfits.open(simgalsets+gsum+str(x)+'.fits')[0].data

            sexinput = image+gimages
            cut3 = glob.glob(simimages+mc+'*.fits')
            if simimages+mc+str(x)+'.fits' in cut3:
                os.remove(simimages+mc+str(x)+'.fits')
            pyfits.writeto(simimages+mc+str(x)+'.fits', sexinput)

            # Run new image (mc) through sextractor
            cmnd = sexdir+'sex '+wdir+simimages+mc+str(x)+'.fits -c '+sexfile
            cmnd += ' -mag_zeropoint '+str(zp)+' -catalog_name'
            cmnd += ' '+wdir+sexcat+sex+str(x)+'.cat'
            os.system(cmnd)

            # Find mag of gsum images (used for aperture correction)
            # using the same apertures from the corresponding mc image
            cmnd2 = sexdir+'sex '+wdir+simimages+mc+str(x)+'.fits,'

            cmnd2 += ''+wdir+simgalsets+gsum+str(x)+'.fits -c defaultsim.sex'
            cmnd2 += ' -mag_zeropoint '+str(zp)+' -catalog_name'
            cmnd2 += ' '+wdir+simcat+sim+str(x)+'.cat'
            os.system(cmnd2)

            if verbose:
                output.write('Generated sextractor catalog: '
                             + sex + str(x) + '.cat \n')

            # Read new sextractor catalog for mc images into table
            sc = Table.read(
                wdir+sexcat+sex+str(x)+'.cat',
                format='ascii.sextractor')

            xpix = sc[sc.colnames[4]]
            ypix = sc[sc.colnames[5]]
            mag = sc[sc.colnames[1]]
            magerr = sc[sc.colnames[2]]
            flux = sc[sc.colnames[3]]
            ellipt = sc[sc.colnames[6]]
            isoarea = sc[sc.colnames[7]]
            flags = sc[sc.colnames[8]]

            # Read in sextractor catalog for gsum images
            simc = Table.read(
                wdir+simcat+sim+str(x)+'.cat',
                format='ascii.sextractor')

            simmag = simc[simc.colnames[1]]
            #simmagerr = simc[simc.colnames[2]]
            #simflux = simc[simc.colnames[3]]

            # Define blank array in which to store sim gal and sex matches
            xmatch = np.zeros(iter[1])
            ymatch = np.zeros(iter[1])
            magmatch = np.zeros(iter[1])
            magerrmatch = np.zeros(iter[1])
            fluxmatch = np.zeros(iter[1])
            elliptmatch = np.zeros(iter[1])
            flagsmatch = np.zeros(iter[1])
            isoareamatch = np.zeros(iter[1])

            simmagmatch = np.zeros(iter[1])
            #simmagerrmatch = np.zeros(iter[1])
            #simfluxmatch = np.zeros(iter[1])

            # Match batch of sim gal with sextractor detections
            for z in xrange(iter[1]):

                # Match x,y elements from simulated with sextractor detections
                xbin = np.asarray(np.where(np.logical_and(
                    xpix >= gxpix[x * iter[1] + z] - pixres,
                    xpix <= gxpix[x * iter[1] + z] + pixres)))[0]

                ybin = np.asarray(np.where(np.logical_and(
                    ypix >= gypix[x * iter[1] + z] - pixres,
                    ypix <= gypix[x * iter[1] + z] + pixres)))[0]

#==============================================================================
#                 # Flag multiple detections
#                 if len(np.intersect1d(xbin,ybin)) == 0:
#                     ymatch[z] += 0
#                     magmatch[z] += 0
#                     magerrmatch[z] += 0
#                     reffmatch[z] += 0
#                     elliptmatch[z] += 0
#                     flagsmatch[z] += 0
#==============================================================================

                # Set x and y positions found by sextractor to variables
                if len(np.intersect1d(xbin, ybin)) == 1:
                    xmatch[z] = xpix[np.intersect1d(xbin, ybin)[0]]
                    ymatch[z] = ypix[np.intersect1d(xbin, ybin)[0]]
                    magmatch[z] = mag[np.intersect1d(xbin, ybin)[0]]
                    magerrmatch[z] = magerr[np.intersect1d(xbin, ybin)[0]]
                    fluxmatch[z] = flux[np.intersect1d(xbin, ybin)[0]]
                    elliptmatch[z] = ellipt[np.intersect1d(xbin, ybin)[0]]
                    flagsmatch[z] = flags[np.intersect1d(xbin, ybin)[0]]
                    isoareamatch[z] = isoarea[np.intersect1d(xbin, ybin)[0]]

                    simmagmatch[z] = simmag[np.intersect1d(xbin, ybin)[0]]
                    #simmagerrmatch[z] = simmagerr[np.intersect1d(xbin,ybin)[0]]
                    #simfluxmatch[z] = simflux[np.intersect1d(xbin,ybin)[0]]
                # Flag multiple detections
                if len(np.intersect1d(xbin, ybin)) > 1:
                    xmatch[z] += -1
                    ymatch[z] += -1
                    magmatch[z] += -1
                    magerrmatch[z] += -1
                    fluxmatch[z] += -1
                    elliptmatch[z] += -1
                    flagsmatch[z] += -1
                    isoareamatch[z] += -1

                    simmagmatch[z] += -1
                    #simmagerrmatch[z] += -1
                    #simfluxmatch[z] += -1

            # Append data to aggregate list
            sxpix.extend(xmatch)
            sypix.extend(ymatch)
            smag.extend(magmatch)
            smagerr.extend(magerrmatch)
            sflux.extend(fluxmatch)
            sellipt.extend(elliptmatch)
            sflags.extend(flagsmatch)
            sisoarea.extend(isoareamatch)

            ssimmag.extend(simmagmatch)
            #ssimmagerr.extend(simmagerrmatch)
            #ssimflux.extend(simfluxmatch)

            # Write data to incremental tables
            f = Table()
            f['simmagmatch'] = simmagmatch
            f['magmatch'] = magmatch
            f['xmatch'] = xmatch
            f['ymatch'] = ymatch
            #f['magerrmatch'] = magerrmatch
            #f['fluxmatch'] = fluxmatch
            #f['elliptmatch'] = elliptmatch
            f['flagsmatch'] = flagsmatch

            f.write(wdir+sexmatch+outfile+str(x)+'.tab', format='ascii.tab')

            ### Insert galfitsersic call here

            t = Table()
            t['ID'] = gid
            t['mag'] = gmag
        #    t['r_eff'] = greff
        #    t['ellipt'] = gellipt
        #    t['pos_angle'] = gpa
        #    t['sxpix'] = list(sxpix)
        #    t['sypix'] = list(sypix)
            t['simmag'] = list(ssimmag)
            t['sexmag'] = list(smag)
        #    t['smag-err'] = list(smagerr)
        #    t['simmag-err'] = list(simmagerr)
        #    t['sr_eff'] = list(sreff)
        #    t['sellipt'] = list(sellipt)
            t['xpix'] = gxpix
            t['ypix'] = gypix
            t['nsersic'] = gsersic
            t['isoarea'] = list(sisoarea)
            t['sflags'] = list(sflags)

            t.write(wdir+'all'+outfile+str(w)+'.tab', format='ascii.tab')
            if verbose:
                output.write('\n data written to all'
                +outfile+str(w)+'.tab \n')


# Vestigial code
'''
    # Setup user options
    usage = "usage: %prog [options] \n"
    usage += "type %prog -h for help"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out",
                      dest="outfilename",
                      default="sexmatch.tab",
                      help="Output to FILE [default : %default]",
                      metavar="FILE")
    parser.add_option("-v", "--verbose",
                      dest="verbose", type=int,
                      default=1,
                      help="Set verbose level to INT [default : %default]",
                      metavar="INT")
    parser.add_option("-i", "--iter", dest="iter", type=int, nargs=2,
                      default=(1000, 10),  # arg1:overall iter; arg2:ngal per image
                      help="Iteration range [default : %default]",
                      metavar="iMin iMax")
    parser.add_option("-x", "--sex", dest="sex",
                      default="default.sex",
                      help="Sextractor file [default : %default]",
                      metavar="FILE")

    (options, args) = parser.parse_args()
'''