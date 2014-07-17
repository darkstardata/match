"""
Created at 4:20 PM on 16 Jul 2014

Project: match
Subprogram: wht2rms

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""
import os
import pyfits                               # Open fits files
import numpy as np

os.chdir('/Users/andrew/iraf')

filter = 'f814w'
inst = 'acs'

# Scan directory for weight maps
image = 'hlsp_candels_hst_'+inst+'_uds-tot_'+filter+'_v1.0_wht.fits'

weight_image = pyfits.open(image, memmap=True)[0].data
rms_image = np.float64(1.0)/np.sqrt(weight_image)
pyfits.writeto('hlsp_candels_hst_'+inst+'_uds-tot_'+filter+'_v1.0_rms.fits', rms_image)
