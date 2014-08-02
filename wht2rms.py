"""
Created at 4:20 PM on 16 Jul 2014

Project: match
Subprogram: wht2rms

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
import os
import sys
import pyfits                               # Open fits files
import numpy as np

os.chdir('/Users/andrew/iraf')
#os.chdir('/home/lokiz/Desktop/fits')


# Field and Filter selector
y = 0       # Field index
z = 1       # Filter index

fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
output_pixelres = [0.06, 0.06, 0.03, 0.03, 0.03]
input_pixelres = [0.13, 0.13, 0.05, 0.05, 0.04]
field = fields[y]
filt = filters[z]

# Instrument selector
if filt is 'f160w' or filt is 'f125w':
    inst = 'wfc3'

elif filt is 'f814w' or filt is 'f606w':
    inst = 'acs'

elif filt is 'f350lp':
    inst = 'uvis'

else:
    sys.exit('Error: No Filter Selected!')


# Scan directory for weight maps
ima = 'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_wht.fits'

weight = pyfits.open(image, memmap=True)[0]
weight_image = weight.data
hdr = weight.header

weight_sqrt = np.sqrt(weight_image)
# Prevents divide by zero
weight_sqrt[weight_sqrt == 0] = 10**(-1*30.0)
rms = np.float32(1.0/weight_sqrt)
print rms[0][0]
rms_image = pyfits.PrimaryHDU(rms, header=hdr)
rms_image.writeto('hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_rms.fits', clobber=True)
