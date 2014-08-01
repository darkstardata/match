"""
Created at 10:35 PM on 05 Jul 2014

Project: match
Subprogram: plotmag

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""

# Graph with title, axes labels, with edge and corner cuts due to noise
import os
import sys
import pyfits
import numpy as np
import collections
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.stats import norm
from astropy.table import Table

wdir = os.getcwd()+'/'

# Field and Filter selector
x = 0       # Field index
y = 0       # Filter index

fields = ['uds', 'bootes']
filters = ['f160w', 'f125w', 'f814w', 'f606w', 'f350lp']
zeropt = [25.96, 26.25, 25.94333, 26.49113, 26.94]
pixel_scale = [0.06, 0.06, 0.03, 0.03, 0.03]


field = fields[x]
filt = filters[y]
zp = zeropt[y]    # zeropoint magnitude
pixscl = pixel_scale[y]    # "/pix

# Instrument selector
if filt is 'f160w' or filt is 'f125w':
    inst = 'wfc3'

elif filt is 'f814w' or filt is 'f606w':
    inst = 'acs'

elif filt is 'f350lp':
    inst = 'uvis'

else:
    sys.exit('Error: No Filter Selected!')

# Local Subdirectories and special input files
science = wdir + 'science/' + field + '/' + filt + '/'

image = pyfits.open(science+'hlsp_candels_hst_'+inst+'_'+field+'-tot_'+filt+'_v1.0_drz.fits')
hdu = image[0]
data = hdu.data
hdr = hdu.header
xsize = hdr['NAXIS1']
ysize = hdr['NAXIS2']


### Cut adjustments
#Cut range for sersic index
nsersicmin = 0
nsersicmax = 12

#Cut range for surface brightness cut (min=19.56, max = 25.24)
sbmin = 18.00
sbmax = 28.00

# Cut galaxies with diameter (pixdia) of galaxies on segmentation map
# Must be an odd number!!!)
pixdia = 3


edge = 10  # edge width around image to exclude outliers
box = 30  # box in top corners to exclude due to heavy noise
#wdir = '/Users/andrew/Documents/Research/Naveen/magiso/'
data = Table.read(science+'allsexmatch.tab', format='ascii.tab')
xid = data[data.colnames[0]]
xpix = data[data.colnames[1]]
ypix = data[data.colnames[2]]
mag = data[data.colnames[3]]
reff = data[data.colnames[4]]
nsersic = data[data.colnames[5]]
sexmag = data[data.colnames[10]]
simmag = data[data.colnames[12]]
npix = data[data.colnames[16]]

magdiffsex = sexmag-mag
magdiffsim = simmag-mag
magdiffaper = sexmag-simmag

sb = np.zeros(np.size(npix))
for w in xrange(np.size(npix)):
    if npix[w] != 0:
        sb[w] = sexmag[w]+2.5*np.log10((pixscl**2)*npix[w])


# Identify overlap between objects inside edge cut and inside corner cuts
xe = np.asarray(np.where(np.logical_and(xpix <= xsize - edge, xpix >= edge)))[0]
ye = np.asarray(np.where(np.logical_and(ypix <= ysize - edge, ypix >= edge)))[0]
elem1 = np.intersect1d(xe, ye)
xp = xpix[elem1]
yp = ypix[elem1]
xc = np.asarray(np.where(np.logical_or(xpix > xsize -  (box + 15), xpix < box)))[0]
yc = np.asarray(np.where(ypix > ysize - box))[0]
elem2 = np.intersect1d(xc, yc)

elem1_multiset = collections.Counter(elem1)
elem2_multiset = collections.Counter(elem2)
overlap = list((elem1_multiset & elem2_multiset).elements())
elem1_remainder = list((elem1_multiset - elem2_multiset).elements())
elem2_remainder = list((elem2_multiset - elem1_multiset).elements())

# uncut portion is elem1_remainder
cutedge = elem1_remainder

# Cut by removing outliers on the edge of the image
binmagdiff = np.asarray(np.where(np.logical_and(magdiffsex >= -2, magdiffsex<=2)))[0]
cutmagdiff = np.intersect1d(binmagdiff, cutedge)

# Find galaxies that fall on other galaxies in segmentation map within
# given radius (pixrad).
# d is the diameter, only work if d is an od number!
# R is radius to search within, r is actual distance form central pixel.
d = pixdia
R = (d-1)/2

# y axis comes first in fits images the elements of segm
segm = pyfits.open(science+'segm/segm_f160w.fits')[0].data
segdet = np.zeros(xpix.size)
for z in xrange(xpix.size):
    xcen = round(xpix[z])
    ycen = round(ypix[z])
    for x in xrange(d):
        for y in xrange(d):
            r = np.sqrt(float(x-R)**2+float(y-R)**2)
            if r < R:
                xelem = int(xcen + x - R)
                if 0 <= xelem <= (xsize-1):
                    yelem = int(ycen + y - R)
                    if 0 <= yelem <= (ysize-1):
                        segdet[z] += segm[yelem, xelem]


# Cut by overlap radius(cut sim galaxies that fall on other galaxies in segmentation map)
binseg = np.asarray(np.where(segdet == 0))[0]
cutseg = np.intersect1d(binseg, cutmagdiff)

# Cut by sersic index
binnsersic = np.asarray(np.where(np.logical_and(nsersic >= nsersicmin, nsersic <= nsersicmax)))[0]
cutnsersic = np.intersect1d(binnsersic, cutseg)

# Cut by surface brightness
binsb = np.asarray(np.where(np.logical_and(sb >= sbmin, sb <= sbmax)))[0]
cutsb = np.intersect1d(binsb, cutnsersic)


# Define precut (before statistics)
prestatcut = cutsb


# best fit of data
(mu, sigma) = norm.fit(magdiffaper[prestatcut])
# 3sigma magnitude cut for graph
binsig = np.asarray(np.where(np.logical_and(magdiffaper >= (mu-5*sigma), magdiffaper <= (mu+5*sigma))))[0]
cutsig = np.intersect1d(binsig, prestatcut)

# Define final cut (after statistics and 3*sigma cuts)
final = cutsig


### Graph Data ###

# Set output image size/quality
fig = plt.figure(dpi=120)

# Define Subplots
gs = gridspec.GridSpec(3, 1, width_ratios=[1], height_ratios=[1, 1, 1], hspace=1.5)
ax1 = plt.subplot(gs[0, :])
ax2 = plt.subplot(gs[-2:, :])

# Plot Histogram of aperture correction
n, bins, patches = ax1.hist(magdiffaper[prestatcut], 100, normed=1)
y = mlab.normpdf(bins, mu, sigma)
ax1.plot(bins, y, 'r--', linewidth=2)
ax1.set_xlabel(r'$\Delta m_{aper} (= m_{sex}-m_{sim})$', fontsize=15)
ax1.set_ylabel(r'Normalized N', fontsize=15)
ax1.set_title(r'MAGISO Aperture Correction: $\mu=%.2f,\ \sigma=%.2f$' % (mu, sigma))
ax1.tick_params(
    axis='x',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',          # ticks along the top edge are off
    labelbottom='on',   # labels along the bottom edge are off
    labelsize=6)     # Size of ticklabels

ax1.tick_params(
    axis='y',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    left='on',          # ticks along the left edge are off
    right='off',        # ticks along the right edge are off
    labelleft='on',     # labels along the left edge are off
    labelsize=6)     # Size of ticklabels

# Plot Scatter of delta mag_aper vs mag_sim
ax2.scatter(simmag[final], magdiffaper[final]-mu, s=5)
ax2.set_xlabel(r'$m_{sim}$', fontsize=15)
ax2.set_ylabel(r'$\Delta m_{aper} (=m_{sex}-m_{sim})$', fontsize=15)
ax2.set_title('MAGISO Sextractor Magnitude Recovery')
ax2.tick_params(
    axis='x',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',          # ticks along the top edge are off
    labelbottom='on',   # labels along the bottom edge are off
    labelsize=10)     # Size of ticklabels

ax2.tick_params(
    axis='y',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    left='on',          # ticks along the left edge are off
    right='off',        # ticks along the right edge are off
    labelleft='on',     # labels along the left edge are off
    labelsize=10)     # Size of ticklabels

# print length(ucm)
print len(final)
plt.show()
plt.savefig('magiso', dpi=120)
