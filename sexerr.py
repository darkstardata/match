"""
Created at 1:25 PM on 31 Jul 2014

Project: match
Subprogram: sexerr

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''
usage: sexerr.py [-h] [--sum] N [N ...]

Process some integers.

positional arguments:
 N           an integer for the accumulator

optional arguments:
 -h, --help  show this help message and exit
 --
 --sum       sum the integers (default: find the max)
'''

import os
import sys
import pyfits
import argparse
import numpy as np
from astropy.table import Table


image = pyfits.open(science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_drz.fits')
hdu = image[0]
data = hdu.data
hdr = hdu.headergaltab
ccdgain = hdr['CCDGAIN']
exptime = hdr['EXPTIME']
inst = hdr['INSTRUME']
filt = hdr['FILTER'].trim()
xpix = hdr['NAXIS1']
ypix = hdr['NAXIS2']







# SExtractor and Galfit directories
# galfit='/home/lokiz/bin/galfit'                                                   # Linux location
galfit = '/usr/local/bin/galfit'                                                    # Mac location
sexdir = '/opt/local/bin/'                                                          # sextractor dir


# SExtractor config files and options
defaultsex = 'default.sex'
defaultsim = 'default.sex'
weight_type = 'MAP_RMS'
weight_image = science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_rms.fits'

#
randgal = 'randgal'                                                                 # galfit config file prefix
random = 'random.tab'                                                               # input list of true random numbers
dustimage = ' '                                                                     # Bad pixel mask for galfit


# Prefixes for iterated input and output files
outfile = 'sexmatch'
galcomb = 'galcomb'                                   # image with n sim gal added
galmodel = 'galmodel'                                 # sim galaxy file prefix
galmodelset = 'galmodelset'
sim = 'sim'
sex = 'sex'

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



# Galgen.py

### Galfit parameters ###



fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
zeropoint = [25.96, 26.25, 25.94333, 26.49113, 26.94]
psf_fwhm = [0.18, 0.12, 0.09, 0.08, 0.08]              # PSF FWHM in units of arcseconds
pixel_scale = [0.06, 0.06, 0.03, 0.03, 0.03]

field = fields[y]
filt = filters[z]
zp = zeropoint[z]
psf = psf_fwhm[z]
pixscl = pixel_scale[z]


psfimage = wdir+'psf/'+'psf_'+filt+'.fits'


# Set outfile name, and verbosity
xpixels_uds = [30720, 30720, 61440, 61440, 61440]
ypixels_uds = [12800, 12800, 25600, 25600, 25600]

xpixels_bootes = 1089
ypixels_bootes = 963

xxpix = u'%7s' % xpix
yypix = u'%7s' % ypix



# Read in sim gal data catalog
galtab = 'galmodel_'+field+'_'+filt+'.tab'
gal = Table.read(science + galtab, format='ascii.tab')
gid = gal[gal.colnames[0]]
gxpix = gal[gal.colnames[1]]
gypix = gal[gal.colnames[2]]
gmag = gal[gal.colnames[3]]
greff = gal[gal.colnames[4]]
gsersic = gal[gal.colnames[5]]
gba = gal[gal.colnames[6]]
gpa = gal[gal.colnames[7]]


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


def add_pstamp(input_image, postage_stamp, xloc, yloc):
    """
    Adds the image 'postage_stamp' to a blank image (zero array) the same size as the 'input_image'.
    The postage stamp size is 300pix X 300pix and as a result the psf is centered on x,y = (151,151).
    Pixel (xcent, ycent) on the 'postage_stamp' is added to pixel (xpos, ypos) in the 'input_image'.
    The remaining pixels are mapped to the surrounding pixels.

    Requires: sys, numpy, pyfits

    Note: FITS images have x and y inverted so y comes first (i.e. [y][x])
    """
    # Find dimensions of 'input_image' and 'postage_stamp'
    ysize_in = input_image.shape[0]
    xsize_in = input_image.shape[1]
    ysize_ps = postage_stamp.shape[0]
    xsize_ps = postage_stamp.shape[1]

    # Fix due to image pixels starting at 1 but array indices starting at 0
    ypos = yloc - 1
    xpos = xloc - 1

    # Define xcent, and ycent by checking if postage stamp size is even or odd number of pixels
    # in x and y dimensions, exit if even since not sure what that will do yet and don't care
    if xsize_ps % 2 == 0:
        xcent = (xsize_ps/2 + 1.0) - 1
    else:
        #xcent = (xsize_ps/2 + 0.5) - 1
        sys.exit(' Error, code not yet ready for postage stamps with odd number of pixels in any dimension!')
    if ysize_ps % 2 == 0:
        ycent = (ysize_ps/2 + 1.0) - 1
    else:
        #ycent = (ysize_ps/2 + 0.5) - 1
        sys.exit(' Error, code not yet ready for postage stamps with odd number of pixels in any dimension!')

    # Make sure 'postage_stamp' dimensions are smaller then the resize image dimensions
    if xsize_ps >= xsize_in or ysize_ps >= ysize_in:
        sys.exit('Error: Postage stamp resizing can only enlarge, ' +
                 'but postage stamp is larger than resized dimension(s)!')
    # Make sure xpos and ypos exist inside resized image dimensions
    if xpos < 0 or xpos > xsize_in-1:
        sys.exit('Error: X center position is outside bounds of input_image dimensions!')
    if ypos < 0 or ypos > ysize_in-1:
        sys.exit('Error: Y center position is outside bounds of input_image dimensions!')

    # Define xmin, xmax, ymin, ymax, xr_offset, xl_offset, yb_offset, and yt_offset
    xr_offset = 0
    xl_offset = 0
    yb_offset = 0
    yt_offset = 0
    ymin = ypos - ycent
    ymax = ypos + ycent - 1
    xmin = xpos - xcent
    xmax = xpos + xcent - 1

    # Check if postage stamp goes outside bounds of image and correct for this.
    if xmin < 0:
        xl_offset = -xmin
        xmin = 0

        if ymin < 0:
            yb_offset = -ymin
            ymin = 0
        elif ymax > ysize_in - 1:
            yt_offset = ymax - (ysize_in - 1)
            ymax = ysize_in - 1

    elif xmax > xsize_in - 1:
        xr_offset = xmax - (xsize_in - 1)
        xmax = xsize_in - 1

        if ymin < 0:
            yb_offset = -ymin
            ymin = 0
        elif ymax > ysize_in - 1:
            yt_offset = ymax - (ysize_in - 1)
            ymax = ysize_in - 1

    else:
        if ymin < 0:
            yb_offset = -ymin
            ymin = 0
        elif ymax > ysize_in - 1:
            yt_offset = ymax - (ysize_in - 1)
            ymax = ysize_in - 1

    ### Add region of postage stamped centered on (xcent, ycent) to
    ### resized image of zero pixels region centered on (xpos, ypos)
    input_image[ymin:ymax+1, xmin:xmax+1] += postage_stamp[yb_offset:ysize_ps-yt_offset, xl_offset:xsize_ps-xr_offset]

    return input_image


# Generate simulated galaxies
def galgen():

    # Output postage stamp dimensions
    xout = 300
    yout = 300

    # Convolution box dimensions
    xconv = 300
    yconv = 300

    # Create lists to store iterated data & sextractor detection data
    txpix = []
    typix = []
    tmag = []
    treff = []
    tsersic = []
    tba = []
    tpa = []

    # Read in table of random numbers, select column of random numbers, and reshape random number array
    rtab = Table.read(random, format='ascii.tab')
    rdat = rtab[rtab.colnames[0]]
    rand = np.reshape(np.array(rdat), (10000, 7))

    # Generate batch of simulated galaxies (total = iteration)
    for x in xrange(iteration):

        r = rand[x]
        rxpix = round(r[0]*xpix, 0)
        rypix = round(r[1]*ypix, 0)
        if rxpix == 0:
            rxpix = 1.00
        if rypix == 0:
            rypix = 1.00
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
                ((xout/2)+1), ((yout/2)+1), rmag, rreff, rsersic, rba, rpa, psf_image=psfimage)
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
    t['ID'].format = '%5.f'
    t['xpix'].format = '%6.f'
    t['ypix'].format = '%6.f'
    t['mag'].format = '%5.2f'
    t['r_eff'].format = '%4.2f'
    t['n_sersic'].format = '%5.2f'
    t['b/a'].format = '%4.2f'
    t['PA'].format = '%5.2f'

    t.write(wdir+'science/'+field+'/'+filt+'/'+galtab, format='ascii.tab')
    if verbose:
        output.write('\n data written to science/'+field+'/'+filt+'/'+galtab+' \n')

    # Exit when finished
    return


def galcomb(real_image):
    # Loop through batches of simulated galaxies
    for w in xrange(iteration[0]):

        # Read in image data and header from science image
        data = real_image
        galset = np.zeros((ypix, xpix))

        # Loop through the individual simulated galaxies
        for x in xrange(iteration[1]):

            # Open simulated galaxy postage stamp, subtract sky, and zero out values less than zero.
            # Then add sku subtracted simulated galaxy to array of zeros.
            gimage = pyfits.open(simgal + galmodel + str(x+w*iteration[1]) + '.fits')[0].data

            data = add_pstamp(data, gimage, gxpix[x+w*iteration[1]], gypix[x+w*iteration[1]])
            galset = add_pstamp(galset, gimage, gxpix[x+w*iteration[1]], gypix[x+w*iteration[1]])

        # Write just simulated galaxy batch to galmodelset*.fits and the simulared + science to galcomb*.fits
        galset_out = pyfits.PrimaryHDU(data=galset)
        galset_out.writeto(simgalset + galmodelset + str(w) + '.fits', clobber=True)
        data_out = pyfits.PrimaryHDU(data=data, header=hdr)
        data_out.writeto(simcomb + galcomb + str(w) + '.fits', clobber=True)



def galsex(iteration, simcomb, galmodelset):
    for x in xrange(iteration):

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


def match():
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


def main():
    parser = argparse.ArgumentParser(description="Characterize sextractor errors for a particular "
                                                 "astronomical fits image given the image and it's"
                                                 " weight map.")
    parser.add_argument('image', metavar='Input Image', type=str, nargs='+',
                       help='input image')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                       const=sum, default=max,
                       help='sum the integers (default: find the max)')

    parser.add_argument('--gen', action='store_true')
    parser.add_argument('--comb', action='store_true')
    parser.add_argument('--sex', action='store_true')
    parser.add_argument('--match', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--verbose', '-v', action='count')

    args = parser.parse_args()

    iteration = [5, 2000]
    output = sys.stdout
    verbose = 1


    # Parameter ranges
    iteration = 10000
    #reff = [0.28, 2.39]                                     # Bootes fields 0.3 -> 2.6 kpc @ z~1.9
    reff = [0.59, 5.10]                                     # UDS WFC3/IR fields 0.3 -> 2.6 kpc @ z~1.9
    #reff = [1.15, 10.20]                                    # UDS WFC/ACS fields 0.3 -> 2.6 kpc @ z~1.9
    ellipt = [0.0, 0.9]
    mag = [18.0, 28.0]
    sersic = [0.01, 12.5]


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



    if gen == True:


    if comb == True:
        galcom(data)
    if sex == True:
        galsex(iteration[0], ,  )

    if match == True:

    return




if __name__ == "__main__":
    main()


