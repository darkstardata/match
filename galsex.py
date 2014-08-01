"""
Created at 4:54 PM on 03 Jul 2014

Project: match
Subprogram: galsex.py

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''
Runs sextractor in single image mode on galcomb*.fits
and runs in dual image mode on galmodelset*.fits to
compare mag inside MAG_ISO aperture in both images.

Runtime ~ 8 hr  10 min
'''

import os
import sys
import pyfits

#Dashboard
iteration = [5, 2000]
output = sys.stdout
verbose = 1
sexdir = '/opt/local/bin/'                  # sextractor dir

# Prefixes for iterated input and output files
galcomb = 'galcomb'                                   # image with n sim gal added
galmodel = 'galmodel'                                 # sim galaxy file prefix
galmodelset = 'galmodelset'
sim = 'sim'
sex = 'sex'

# Field and Filter selector
y = 0       # Field index
z = 1       # Filter index

fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
zeropoint = [25.96, 26.25, 25.94333, 26.49113, 26.94]
psf_fwhm = [0.18, 0.12, 0.09, 0.08, 0.08]              # PSF FWHM in units of arcseconds
field = fields[y]
filt = filters[z]
zp = zeropoint[z]
psf = psf_fwhm[z]


# Local Subdirectories and special input files
wdir = os.getcwd()+'/'
science = wdir+'science/'+field+'/'+filt+'/'
simcomb = science+'simcomb/'
simgal = science+'simgal/'
sexconfig = science+'sexconfig/'
simgalset = science+'simgalset/'
simcat = science+'simcat/'
sexcat = science+'sexcat/'

galtab = 'galmodel_'+field+'_'+filt+'.tab'


# Instrument selector
if filt is 'f160w' or filt is 'f125w':
    inst = 'wfc3'

elif filt is 'f814w' or filt is 'f606w':
    inst = 'acs'

elif filt is 'f350lp':
    inst = 'uvis'

else:
    sys.exit('Error: No Filter Selected!')


# Sextractor config file parameters
defaultsex = 'default.sex'
defaultsim = 'default.sex'
weight_type = 'MAP_RMS'
weight_image = science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_rms.fits'


# Scan directory for input fits files while excluding program files and alert if no .fits input images are present
# Read in image data and header from science image
image = pyfits.open(science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_drz.fits')
hdu = image[0]
data = hdu.data
hdr = hdu.header
ccdgain = hdr['CCDGAIN']
exptime = hdr['EXPTIME']


for x in xrange(iteration[0]):

    # Run sextractor in single image mode on new image (galcomb*.fits)
    cmnd = sexdir + 'sex ' + simcomb + galcomb + str(x) + '.fits'
    cmnd += ' -c ' + sexconfig + defaultsex
    cmnd += ' -catalog_name ' + sexcat + sex + str(x) + '.cat'
    cmnd += ' -GAIN ' + str(ccdgain * exptime)
    cmnd += ' -MAG_ZEROPOINT ' + str(zp)
    cmnd += ' -SEEING_FWHM ' + str(psf)
    cmnd += ' -WEIGHT_TYPE ' + weight_type
    cmnd += ' -WEIGHT_IMAGE ' + weight_image
    os.system(cmnd)

    if verbose:
        output.write('Generated sextractor catalog: ' + sex+str(x)+'.cat \n')

    # Run sextractor in dual image mode:
    #   Find mag of galmodelset*.fits image (used for aperture correction)
    #   using the same apertures from the corresponding galcomb*.fits image
    cmnd2 = sexdir + 'sex ' + simcomb + galcomb + str(x) + '.fits,'
    cmnd2 += simgalset + galmodelset + str(x) + '.fits'
    cmnd2 += ' -c ' + sexconfig + defaultsex
    cmnd2 += ' -catalog_name ' + simcat + sim + str(x) + '.cat'
    cmnd2 += ' -GAIN ' + str(ccdgain * exptime)
    cmnd2 += ' -MAG_ZEROPOINT ' + str(zp)
    cmnd2 += ' -SEEING_FWHM ' + str(psf)
    cmnd2 += ' -WEIGHT_TYPE ' + weight_type +',NONE'
    cmnd2 += ' -WEIGHT_IMAGE ' + weight_image +','+ weight_image
    cmnd2 += ' -WEIGHT_GAIN N,N'
    cmnd2 += ' -BACK_TYPE MANUAL,MANUAL'
    cmnd2 += ' -BACK_VALUE 0.0,0.0'
    os.system(cmnd2)

    if verbose:
        output.write('Generated sextractor catalog: ' + sim+str(x)+'.cat \n')