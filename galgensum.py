# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 12:25:30 2014
Fixed on Mon Jul  2 13:19:56 2014

@author: andrew
"""
'''
Combines single simulated galaxy images into batch images with the number of galaxies equal to iteration[1]

Runtime ~ 8 minutes
'''
import os
import sys
import glob                                 # Find files
import numpy as np
import pyfits                               # Open fits files


#Dashboard
xx = 1089
yy = 963


# Subdirectories
simgal = 'simgal/'              # Dir for individual Galaxies
simgalset = 'simgalset/'      # Dir for set of galaxies to be added to image


# Prefixes for output files
galmodel = 'galmodel'           # Prefix for individual Galaxies
galmodelset = 'galmodelset'                   # Prefix for set of galaxies to be added to image
mc = 'mc'                       # Prefix for set of galaxies plus image


iteration = [1000, 10]          # [ number of output images, number of galaxies per image]
verbose = 1
output = sys.stdout


def main():

    # Loop through iteration[0] batches
    for x in xrange(iteration[0]):

        # Define blank array for image
        gimages = np.zeros((yy, xx))

        # Check dir for file with same name as output, remove if present
        cut = glob.glob(simgalset+galmodelset+'*.fits')
        if simgalset+galmodelset+str(x)+'.fits' in cut:
            os.remove(simgalset+galmodelset+str(x)+'.fits')
            
        # Loop through iteration[1] per batch
        for y in xrange(iteration[1]):
                
            # Add single simulated galaxy to the batch image
            gimage = pyfits.open(simgal+galmodel+str(x*iteration[1]+y)+'.fits')
            gimages += gimage[0].data

            # Print message when each image merger is complete
            if x == (iteration[1]-1):
                if verbose is 1:
                    output.write('\n ' + str(iteration[1]) + ' simulated galaxies merged into a ' +
                                 'single image; ' + galmodelset + str(x * iteration[1] + y) + '.fits \n')

        # Check dir for file with same name as output, remove if present
        cut = glob.glob(simgalset+galmodelset+'*.fits')
        if simgalset+galmodelset+str(x)+'.fits' in cut:
            os.remove(simgalset+galmodelset+str(x)+'.fits')

        # Write image (simulated batch plus science) to galmodelset*.fits (* is the iteration number)
        pyfits.writeto(simgalset+galmodelset+str(x)+'.fits', gimages)

        # Print message when each image merger is complete
        if x == (iteration[1]-1):
            if verbose is 1:
                output.write('\n '+str(iteration[1])+' simulated galaxies merged into a single image;'+
                             ' galmodelset'+str(x*iteration[1]+y)+'.fits \n')

if __name__ == "__main__":
   main()