# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 12:25:30 2014

@author: andrew
"""
import os
import glob                                 # Find files
import numpy as np
import pyfits                               # Open fits files
from optparse import OptionParser           # Parse options

#Dashboard
xx = 1089
yy = 963

# Set working directory
wdir = os.getcwd()+'/'

# Subdirectories
simgal = 'simgal/'              # Dir for individual Galaxies
simgalsets = 'simgalsets/'      # Dir for set of galaxies to be added to image
simimages = 'simimages/'        # Dir for set of galaxies plus image

# Prefixes for output files
galmodel = 'galmodel'           # Prefix for individual Galaxies
gsum = 'gsum'                   # Prefix for set of galaxies to be added to image
mc = 'mc'                       # Prefix for set of galaxies plus image

def main():
    # Setup user options
    usage = "usage: %prog [options] \n"
    usage+= "type %prog -h for help"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", 
                      dest="verbose", type=int, 
                      default=1,
                      help="Set verbose level to INT [default : %default]",
                      metavar="INT")
    parser.add_option("-i", "--iter", dest="iter", type=int, nargs=2,
                      default=(1000, 10), #arg1:overall iter; arg2:ngal per image
                      help="Iteration range [default : %default]",
                      metavar="iMin iMax")

    
    (options, args) = parser.parse_args()

    # Set outfile name, number of iterations, and verbosity
    iter=options.iter
 

    for x in xrange(iter[0]):
        gimages = np.zeros((yy,xx))
        
        # Loop through the ngal(i.e.iter[1]) per batch to add to image
        for y in xrange(iter[1]):
                
            # Load single sim gal image
            gimage = pyfits.open(
                simgal+galmodel+str(x*iter[1]+y)+'.fits')
            
            # Combine sim gals into batch sim gal image
            gimages += gimage[0].data                   
        
        
        # Write input fits image to run through sextractor

            
        cut = glob.glob(simgalsets+gsum+'*.fits')
        if simgalsets+gsum+str(x)+'.fits' in cut: 
            os.remove(simgalsets+gsum+str(x)+'.fits')
            
        pyfits.writeto(simgalsets+gsum+str(x)+'.fits',gimages)
        
                
if __name__ == "__main__":
   main()