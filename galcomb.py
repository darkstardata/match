"""
Created at 7:56 PM on 02 Jul 2014

Project: match
Subprogram: galcomb

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
'''
Combines science images with galmodel*.fits images and outputs galcomb*.fits for processing with sextractor.
Also outputs an image with just the simulated galaxies that were added in each batch called simgalset*.fits .

Runtime ~ 1 hr  5 min
'''
import os
import sys
import pyfits                               # Open fits files
import numpy as np
from astropy.table import Table             # Open simulated galaxy data table


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


def main():
    #Dashboard
    wdir = os.getcwd()+'/'

    # Prefixes for iterated input and output files
    galcomb = 'galcomb'                                   # image with n sim gal added
    galmodel = 'galmodel'                                 # sim galaxy file prefix
    galmodelset = 'galmodelset'

    # Field and Filter selector
    y = 0       # Field index
    z = 2       # Filter index

    fields = ['uds', 'bootes']
    filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']

    field = fields[y]
    filt = filters[z]

    # Local Subdirectories and special input files
    science = wdir + 'science/' + field + '/' + filt + '/'
    simcomb = science + 'simcomb/'
    simgal = science + 'simgal/'
    simgalset = science + 'simgalset/'

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

    # Read and array data from simulated galaxy table (galmodel.tab)
    gal = Table.read(science + galtab, format='ascii.tab')
    gxpix = gal[gal.colnames[1]]
    gypix = gal[gal.colnames[2]]

    # Set iterations and name input image
    iteration = [5, 2000]                  # [ number of combined images, number of sim galaxies per combined image]

    # Loop through batches of simulated galaxies
    for w in xrange(iteration[0]):

        # Read in image data and header from science image
        image = pyfits.open(science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_drz.fits')
        hdu = image[0]
        data = hdu.data
        hdr = hdu.header
        xpix = hdr['NAXIS1']
        ypix = hdr['NAXIS2']
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


if __name__ == "__main__":
    main()
