"""
Created on Sat Jan  4 19:42:20 2014
Fixed on Mon Jun 30 21:15:35 2014

Project: match
Subprogram: galgen

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside
"""

"""
This program generates a batch of simulated galaxies written to fits files
as well as a catalog file of the simulated galaxies properties which include 
7 parameters;

    x position (pixels)
    y position (pixels)
    magnitude
    effective radius (half-light radius)
    ellipticity
    sersic index (n)
    position angle.

Assumed Cosmology: h = 0.71, Omega_M = 0.27, Omega_vac = 0.73, z = 1.9.
(for determining proper range for effective radii in pixels using pixel scale and cosmology)

Runtime ~ 14 hr
"""

import os
import sys
import numpy as np
from astropy.table import Table


#Dashboard

wdir = os.getcwd()+'/'
# galfit='/home/lokiz/bin/galfit'                                                   # Linux location
galfit = '/usr/local/bin/galfit'                                                    # Mac location
randgal = 'randgal'                                                                 # galfit config file prefix
galmodel = 'galmodel'                                                               # sim galaxy file prefix
random = 'random.tab'                                                               # input list of true random numbers
dustimage = ' '

# Set outfile name, and verbosity
outfilename = 'galmodel.tab'
verbose = 1
output = sys.stdout

# Output postage stamp dimensions
xout = 300
yout = 300

# Convolution box dimensions
xconv = 300
yconv = 300

# Parameter ranges
iteration = 10000
#reff = [0.28, 2.39]                                     # Bootes fields 0.3 -> 2.6 kpc @ z~1.9
reff = [0.59, 5.10]                                     # UDS WFC3/IR fields 0.3 -> 2.6 kpc @ z~1.9
#reff = [1.15, 10.20]                                    # UDS WFC/ACS fields 0.3 -> 2.6 kpc @ z~1.9
ellipt = [0.0, 0.9]
mag = [18.0, 28.0]
sersic = [0.01, 12.5]


# Field and Filter selector
w = 0       # Field index
z = 0       # Filter index

fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
zeropt = [25.96, 26.25, 25.94333, 26.49113, 26.94]
pixel_scale = [0.06, 0.06, 0.03, 0.03, 0.03]

field = fields[w]
filt = filters[z]
zp = zeropt[z]    # zeropoint magnitude
pixscl = pixel_scale[z]    # "/pix
psfimage = wdir+'psf/'+'psf_'+filt+'.fits'

xpixels_uds = [30720, 30720, 61440, 61440, 61440]
ypixels_uds = [12800, 12800, 25600, 25600, 25600]

xpixels_bootes = 1089
ypixels_bootes = 963

if field is 'uds':
    xpix = xpixels_uds[z]
    ypix = ypixels_uds[z]

elif field is 'bootes':
    xpix = xpixels_bootes
    ypix = ypixels_bootes

else:
    sys.exit('Error: No Field Selected!')

xxpix = u'%7s' % xpix
yypix = u'%7s' % ypix


# Prevent code breaks
if isinstance(iteration, int) is False:
    sys.exit('Interations must be an integer value')
if iteration < 1:
    sys.exit('Minimum interations must be equal to or greater than one!')
if reff[0] <= 0:
    sys.exit('Minimum effective radius must be greater than zero!')
if ellipt[0] < 0:
    sys.exit('Minimum ellipticity must be equal to or greater than zero!')
if ellipt[1] >= 1:
    sys.exit('Maximum ellipticity must be less than one!')
if sersic[0] <= 0:
    sys.exit('Minimum sersic index must be greater than zero!')


# Generate galfit feed file
def galfeed(output_image, xpix_image, ypix_image, xpix_conv, ypix_conv, zeropoint, pixelscale, xpix_coord, ypix_coord,
            magnitude, radeff, n_sersic, ba, pa, psf_image='', dust_image=''):

    # Generate string versions for GALFIT config file. Rounds all parameters to 2 decimal places
    sxpix = u'%7s' % u'%7.2f' % xpix_coord
    sypix = u'%7s' % u'%7.2f' % ypix_coord
    smag = u'%5s' % u'%5.2f' % magnitude
    sreff = u'%8s' % u'%3.2f' % radeff
    ssersic = u'%5s' % u'%5.2f' % n_sersic
    spa = u'%5s' % u'%5.2f' % pa
    sab = u'%7s' % u'%7.2f' % ba   # b/a = 1 - e

    # Creating GALFIT feed file for simulated galaxy
    galfitfeed = open(randgal + '.feed', 'w')
    galfitfeed.write('# IMAGE PARAMETERS\n')
    galfitfeed.write('A) '                    '# Input Data image (FITS file)\n')
    galfitfeed.write('B) '+output_image+'      # Name for the output image\n')
    galfitfeed.write('C) none                  # Noise image name (made from data if blank or "none")\n')
    galfitfeed.write('D) '+psf_image+'          # Input PSF image and (optional) diffusion kernel\n')
    galfitfeed.write('E) 1                     # PSF oversampling factor relative to data\n')
    galfitfeed.write('F) '+dust_image+'         # Pixel mask (ASCII file or FITS file with non-0 values)\n')
    galfitfeed.write('G) none                  # Parameter constraint file (ASCII)\n')
    galfitfeed.write('H) 1 '+str(xpix_image)+' 1 '+str(ypix_image)+'   # Image region to fit (xmin xmax ymin ymax)\n')
    galfitfeed.write('I) '+str(xpix_conv)+'   '+str(ypix_conv)+'             # Size of convolution box (x y)\n')
    galfitfeed.write('J) ' + str(zeropoint) + '       # Magnitude photometric zeropoint\n')
    galfitfeed.write('K) '+str(pixelscale)+' '+str(pixelscale)+'           # Plate scale (dx dy)\n')
    galfitfeed.write('O) regular               # Display type (regular, curses, both)\n')
    galfitfeed.write('P) 0                     # Create ouput only? (1=yes; 0=optimize)\n')
    galfitfeed.write('S) 0                     # Modify/create objects interactively?\n')
    galfitfeed.write('\n')
    galfitfeed.write('\n')
    galfitfeed.write('# Sersic function\n')
    galfitfeed.write('\n')
    galfitfeed.write(' 0) sersic               # Object type\n')
    galfitfeed.write(' 1) ' + sxpix + ' ' + sypix + '  0 0 # position x, y       [pixels]\n')
    galfitfeed.write(' 3) ' + smag + '          0   # total magnitude\n')
    galfitfeed.write(' 4) ' + sreff + '          0   # R_e                 [pixels]\n')
    galfitfeed.write(' 5) ' + ssersic + '          0   # Sersic exponent (deVauc=4, expdisk=1)\n')
    galfitfeed.write(' 6)                  0   # --------\n')
    galfitfeed.write(' 7)                  0   # --------\n')
    galfitfeed.write(' 8)                  0   # --------\n')
    galfitfeed.write(' 9) ' + sab + '          0   # axis ratio (b/a)\n')
    galfitfeed.write('10) ' + spa + '          0   # position angle (PA)  [Degrees: Up=0, Left=90]\n')
    galfitfeed.write(' Z)                  0   # output image (see above)\n')
    galfitfeed.write('\n')
    galfitfeed.close()


# Generate simulated galaxies
def main():

    # Create lists to store iterated data & sextractor detection data
    txpix = []
    typix = []
    tmag = []
    treff = []
    tsersic = []
    tba = []
    tpa = []

    # Read in table of random numbers, select column of random numbers, and reshape random numer array
    rtab = Table.read(random, format='ascii.tab')
    rdat = rtab[rtab.colnames[0]]
    rand = np.reshape(np.array(rdat), (10000, 7))

    # Generate batch of simulated galaxies (total = iteration)
    for x in xrange(iteration):
#   for x in range(8749,9999,1):

        r = rand[x]
        rxpix = round(r[0]*xpix, 2)
        rypix = round(r[1]*ypix, 2)
        rmag = round(r[2]*abs(mag[1]-mag[0])+mag[0], 2)
        rreff = round(r[3]*abs(reff[1]-reff[0])+reff[0], 2)
        rsersic = round(r[5]*abs(sersic[1]-sersic[0])+sersic[0], 2)
        rba = round(1-r[4]*abs(ellipt[1]-ellipt[0])+ellipt[0], 2)
        rpa = round(r[6]*90, 2)

        # Store random parameters in lists that will be added to table t
        txpix.extend([rxpix])
        typix.extend([rypix])
        tmag.extend([rmag])
        treff.extend([rreff])
        tsersic.extend([rsersic])
        tba.extend([rba])
        tpa.extend([rpa])

        # Generate galfit feed file and run galfit
        galfeed(wdir+'science/'+field+'/'+filt+'/simgal/'+galmodel+str(x)+'.fits', xout, yout, xconv, yconv, zp, pixscl,
                xout/2, yout/2, rmag, rreff, rsersic, rba, rpa, psf_image=psfimage)
        cmnd1 = galfit + ' -noskyest ' + randgal + '.feed'
        os.system(cmnd1)

        # Print message when each 1/10th of total iterations are complete
        if isinstance(float(x)*10/iteration, int) is True:
            if verbose is 1:
                output.write('\n '+str(iteration/10)+' simulated galaxies created \n')

    t = Table()
    t['ID'] = list(np.arange(iteration))
    t['xpix'] = txpix
    t['ypix'] = typix
    t['mag'] = tmag
    t['r_eff'] = treff
    t['n_sersic'] = tsersic
    t['b/a'] = tba
    t['PA'] = tpa

    # Format table data
    t['ID'].format = '%4d'
    t['xpix'].format = '%8.2f'
    t['ypix'].format = '%8.2f'
    t['mag'].format = '%5.2f'
    t['r_eff'].format = '%4.2f'
    t['n_sersic'].format = '%5.2f'
    t['b/a'].format = '%4.2f'
    t['PA'].format = '%5.2f'

    t.write(wdir+'science/'+field+'/'+filt+'/'+outfilename, format='ascii.tab')
    if verbose:
        output.write('\n data written to science/'+field+'/'+filt+'/'+outfilename+' \n')


if __name__ == "__main__":
    main()