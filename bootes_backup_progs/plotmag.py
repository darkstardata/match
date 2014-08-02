"""
Created at 10:35 PM on 05 Jul 2014

Project: match
Subprogram: plotmag

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""

# Graph with title, axes labels, with edge and corner cuts due to noise
import numpy as np
import collections
from astropy.table import Table
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import norm
import pyfits
import matplotlib.gridspec as gridspec

n=0


### Cut adjustments
#Cut range for sersic index
nsersicmin = 0
nsersicmax = 2

#Cut range for surface brightness cut (min=19.56, max = 25.24)
sbmin = 19.00
sbmax = 23.00

# Cut galaxies with diameter (pixdia) of galaxies on segmentation map
# Must be an odd number!!!)
pixdia = 3


edge = 10  # edge width around image to exclude outliers
box = 30  # box in top corners to exclude due to heavy noise
#wdir = '/Users/andrew/Documents/Research/Naveen/magiso/'
wdir = '/Users/andrew/git/match/'
data = Table.read(wdir+'allsexmatch'+str(n)+'.tab', format='ascii.tab')
xid = data[data.colnames[0]]
xpix = data[data.colnames[1]]
ypix = data[data.colnames[2]]
mag = data[data.colnames[3]]
nsersic = data[data.colnames[5]]
sexmag = data[data.colnames[10]]
simmag = data[data.colnames[12]]
npix = data[data.colnames[16]]

magdiffsex = sexmag-mag
magdiffsim = simmag-mag

sb = np.zeros(np.size(npix))
for w in xrange(np.size(npix)):
    if npix[w] != 0:
        sb[w] = sexmag[w]+2.5*np.log10((0.128**2)*npix[w])


# Identify overlap between objects inside edge cut and inside corner cuts
xe = np.asarray(np.where(np.logical_and(xpix<=1089-edge, xpix>=edge)))[0]
ye = np.asarray(np.where(np.logical_and(ypix<=963-edge, ypix>=edge)))[0]
elem1 = np.intersect1d(xe,ye)
xp = xpix[elem1]
yp = ypix[elem1]
xc = np.asarray(np.where(np.logical_or(xpix>1089-(box+15), xpix<box)))[0]
yc = np.asarray(np.where(ypix>963-box))[0]
elem2 = np.intersect1d(xc,yc)

elem1_multiset = collections.Counter(elem1)
elem2_multiset = collections.Counter(elem2)
overlap = list((elem1_multiset & elem2_multiset).elements())
elem1_remainder = list((elem1_multiset - elem2_multiset).elements())
elem2_remainder = list((elem2_multiset - elem1_multiset).elements())

# uncut portion is elem1_remainder
cutedge = elem1_remainder

# Cut by removing outliers on the edge of the image
binmagdiff = np.asarray(np.where(np.logical_and(magdiffsex>=-2,magdiffsex<=2)))[0]
cutmagdiff = np.intersect1d(binmagdiff,cutedge)

# Find galaxies that fall on other galaxies in segmentation map within
# given radius (pixrad).
d=pixdia
l=(d-1)/2
segm = pyfits.open(wdir+'segm/segm'+str(n+1)+'0.fits')[0].data
segdet = np.zeros(xpix.size)
for z in xrange(xpix.size):
    xcen = round(xpix[z])
    ycen = round(ypix[z])
    for x in xrange(d):
        for y in xrange(d):
            len = np.sqrt(float(x-l)**2+float(y-l)**2)
            if len<l:
                xelem = int(xcen + x - l)
                if xelem>=0 and xelem<=962:
                    yelem = int(ycen + y - l)
                    if yelem>=0 and yelem<=1088:
                        segdet[z] += segm[xelem,yelem]


# Cut by overlap radius(cut sim galaxies that fall on other galaxies in segmentation map)
binseg = np.asarray(np.where(segdet==0))[0]
cutseg = np.intersect1d(binseg,cutmagdiff)

# Cut by sersic index
binnsersic = np.asarray(np.where(np.logical_and(
            nsersic>=nsersicmin,nsersic<=nsersicmax)))[0]
cutnsersic = np.intersect1d(binnsersic,cutseg)

# Cut by surface brightness
binsb = np.asarray(np.where(np.logical_and(
            sb>=sbmin,sb<=sbmax)))[0]
cutsb = np.intersect1d(binsb,cutnsersic)


# Define precut (before statistics)
prestatcut = cutsb


# best fit of data
(mu, sigma) = norm.fit(magdiffsex[prestatcut])
# 3sigma magnitude cut for graph
binsig = np.asarray(np.where(np.logical_and(magdiffsex>=(
        mu-3*sigma),magdiffsex<=(mu+3*sigma))))[0]
cutsig = np.intersect1d(binsig,prestatcut)

# Define final cut (after statistics and 3*sigma cuts)
final = cutsig


# Graph title
#fig.suptitle(r'2-D posterior probability density - $prob( \alpha , \beta  | \{ x_k \} )$', fontsize=18)

#==============================================================================
# # Subplots
# gs1 = gridspec.GridSpec(3, 3)
# gs1.update(left=0.05, right=0.48, wspace=0.05)
# ax1 = plt.subplot(gs1[:, :])
# gs2 = gridspec.GridSpec(3, 3)
# gs2.update(left=0.55, right=0.98, hspace=0.05)
# ax2 = plt.subplot(gs2[:, :])
#==============================================================================

# Set output image size/quality
fig = plt.figure(dpi=120)

# Define Subplots
gs = gridspec.GridSpec(3,1, width_ratios=[1], height_ratios=[1,1,1], hspace=0.9)
ax1 = plt.subplot(gs[0, :])
ax2 = plt.subplot(gs[-2:, :])

# Plot Histogram of aperture correction
n, bins, patches = ax1.hist(magdiffsex[prestatcut], 100, normed=1)
y = mlab.normpdf( bins, mu, sigma)
ax1.plot(bins, y, 'r--', linewidth=2)
ax1.set_xlabel(r'$\Delta m_{sim} (= simmag-mag)$', fontsize=15)
ax1.set_ylabel(r'Normalized N', fontsize=15)
ax1.set_title(r'MAGISO Aperture Correction: $\mu=%.2f,\ \sigma=%.2f$' %(mu, sigma))
ax1.tick_params(\
    axis='x',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',          # ticks along the top edge are off
    labelbottom='on',   # labels along the bottom edge are off
    labelsize = 6)     # Size of ticklabels

ax1.tick_params(\
    axis='y',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    left='on',          # ticks along the left edge are off
    right='off',        # ticks along the right edge are off
    labelleft='on',     # labels along the left edge are off
    labelsize = 6)     # Size of ticklabels

# Plot Scatter of delta_mag vs mag
ax2.scatter(mag[final], magdiffsex[final]-mu, s=5)
ax2.set_xlabel(r'$m$', fontsize=20)
ax2.set_ylabel(r'$\Delta m_{sex} (=sexmag-mag)$', fontsize=20)
ax2.set_title('MAGISO Sextractor Magnitude Recovery')
ax2.tick_params(\
    axis='x',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',          # ticks along the top edge are off
    labelbottom='on',   # labels along the bottom edge are off
    labelsize = 10)     # Size of ticklabels

ax2.tick_params(\
    axis='y',           # changes apply to the y-axis
    which='both',       # both major and minor ticks are affected
    left='on',          # ticks along the left edge are off
    right='off',        # ticks along the right edge are off
    labelleft='on',     # labels along the left edge are off
    labelsize = 10)     # Size of ticklabels

# print len(ucm)
plt.show()
fig.savefig('magiso_'+str(sbmin)+'_'+str(sbmax)+'', dpi=120)
