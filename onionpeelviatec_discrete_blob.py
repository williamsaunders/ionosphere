# imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import ticker, cm
import matplotlib.patches as mpatches
import copy
from numpy import matrix
import sys

plt.rcParams['font.size']=12
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

rpkm = 3400. # Rmars in km
r0km = 120. + rpkm # altitude of electron density peak at subsolar point in km in radial distance
n0cm3 = 2e5 # peak electron density at the subsolar point, cm-3
dsh0km = 10. # density scale height, km

# create grid to keep track of electron density
nx = 2001
nc = round((nx-1.)/2)
ny = 1001

xarr1dkm = np.linspace(-rpkm, rpkm, nx)
yarr1dkm = np.linspace(rpkm, 2*rpkm, ny)

xarr2dkm = np.zeros((nx, ny))
yarr2dkm = np.zeros((nx, ny))

# populate 2d arrays from 1d arrays 
for i, x in zip(range(nx), xarr1dkm):
    xarr2dkm[i,:] = x
    yarr2dkm[i,:] = yarr1dkm

# transform into polar coodinates
rarr2dkm = np.sqrt(xarr2dkm**2 + yarr2dkm**2)
theta2rad = np.arctan2(yarr2dkm, xarr2dkm)
theta2ddeg = theta2rad * 180./np.pi

# define spherical ionosphere                                                                                                         
dummyx2d = (rarr2dkm - r0km) / dsh0km
nelec2dcm3 = n0cm3 * np.exp(1 - dummyx2d - np.exp(-dummyx2d))

# CURVED BLOB
y300km = np.where(np.isclose(yarr1dkm-rpkm, 300, atol=.5*(yarr1dkm[1]-yarr1dkm[0])))[0][0] #y index of 300 km
mu = (int(nx/2.),y300km) ### CENTERED BLOB
peak = nelec2dcm3[mu[0],mu[1]]

rarr2dkm0 = np.abs(rarr2dkm - (rpkm + 300))
rstretchfactor = 4
rsomething1 = np.exp(-rarr2dkm0*rstretchfactor**-1)
rsomething2 = rsomething1/np.max(rsomething1)
rsomething3 = copy.copy(rsomething2)

tstretchfactor = 1e10
tsomething1 = np.abs(90-theta2ddeg)
tsomething2 = np.exp(-tsomething1**2*tstretchfactor**-1)
tsomething3 = tsomething2/np.max(tsomething1)

nelec2dcm3_g = n0cm3*rsomething3*tsomething3*peak
nelec2dcm3_orig = copy.copy(nelec2dcm3)
nelec2dcm3 = nelec2dcm3 + nelec2dcm3_g


# convert to m-3 
nelec2dm3 = nelec2dcm3 * 1e6
nelec2dm3_g = nelec2dcm3_g * 1e6
nelec2dm3_orig = nelec2dcm3_orig*1e6
n0m3 = n0cm3 * 1e6


# plot just the blob
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log10(nelec2dcm3_g)
x0, y0 = np.where(Z==-np.inf)
Z2 = copy.copy(Z)
Z2[x0,y0] = -70000
contours = [-70000] + np.arange(-30,30,5).tolist()
CS = plt.contourf(X, Y, Z2,  contours, cmap=cm.PuBu_r, vmin=-30, vmax=20)
cbar = plt.colorbar(CS,ticks=contours)
cbar.ax.set_ylabel(r'log electron density [cm$^{-3}$]')
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
plt.tight_layout()
plt.title('blob')
#plt.savefig('/Users/saunders/Documents/planet_research/blob_discrete/blob.png', bbox_inches='tight', dpi=200)
plt.show()


# make contour plot of electorn density
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log10(nelec2dm3)
x0, y0 = np.where(Z==-np.inf)
Z2 = copy.copy(Z)
Z2[x0,y0] = -70000
CS = plt.contourf(X, Y, Z2,  contours, cmap=cm.PuBu_r, vmin=-30, vmax=20)
cbar = plt.colorbar(CS,ticks=contours)
cbar.ax.set_ylabel(r'log electron density [m$^{-3}$]')

# formatting and labeling
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
plt.tight_layout()
#plt.ylim(3400,4400)
#plt.xlim(-500,500)
plt.title('blob+symmetrical')
#plt.savefig('/Users/saunders/Documents/planet_research/blob_discrete/contour_blob.png', bbox_inches='tight', dpi=200)
plt.show()

# unit conversions
xarr1dm = xarr1dkm * 1e3
yarr1dm = yarr1dkm * 1e3

#  ABEL TRANSFORM SOLUTION
dl = xarr1dm[1] - xarr1dm[0]
tec1dm2 = np.sum(nelec2dm3, axis=0) * dl
tec1dm2_orig = np.sum(nelec2dm3_orig, axis=0) * dl

# change variable names 
newbigx = np.array(tec1dm2)
newbigx_orig = np.array(tec1dm2_orig)
xa = np.array(yarr1dm)

# perform derivative: dOmega/dX
newdxdr = (np.roll(newbigx,-1) - np.roll(newbigx,1)) / (np.roll(xa,-1) - np.roll(xa,1))
newdxdr[0] = 2.*newdxdr[2] - newdxdr[1]
newdxdr[-1] = 2.*newdxdr[-3] - newdxdr[-2]

newdxdr_orig  = (np.roll(newbigx_orig,-1) - np.roll(newbigx_orig,1)) / (np.roll(xa,-1) - np.roll(xa,1))
newdxdr_orig[0] = 2.*newdxdr_orig[2] - newdxdr_orig[1]
newdxdr_orig[-1] = 2.*newdxdr_orig[-3] - newdxdr_orig[-2]

# some math checks 
pp = np.argsort(xa) # make sure order is correct becuase sometimes radio occultations \ 
                    # come out of order becuase of refraction or direction of transit 
reva2 = xa[pp] # ordered impact parameter, m (closest approach distance)
thisrkm = reva2 / 1e3 # ordered impact parameter, km
newrevdxdr = newdxdr[pp] # ordered spatial derivative of bigx, m-3
newrevdxdr_orig = newdxdr_orig[pp] 

# set first element of newrevdxdr to 0 to avoid problems with numerical integration
newrevdxdr[-1] = 0.
newrevdxdr_orig[-1] = 0.

# perform integration

invnelec1dm3 = []
invnelec1dm3_orig = []

print('\n')
i = 0
while i < len(newrevdxdr):
    print('\r', i, '/', ny, end='   ')
    stuff = 0.
    stuff_orig = 0.
    j = i
    while j < len(newrevdxdr) - 1:
        val = .5 * ( np.log(reva2[j]/reva2[i] + np.sqrt((reva2[j]/reva2[i])**2. - 1.)) + \
                np.log(reva2[j+1]/reva2[i] + np.sqrt((reva2[j+1]/reva2[i])**2. - 1.)) ) \
                * (newrevdxdr[j+1] - newrevdxdr[j])
        val_orig = .5 * ( np.log(reva2[j]/reva2[i] + np.sqrt((reva2[j]/reva2[i])**2. - 1.)) + \
                     np.log(reva2[j+1]/reva2[i] + np.sqrt((reva2[j+1]/reva2[i])**2. - 1.)) ) \
                     * (newrevdxdr_orig[j+1] - newrevdxdr_orig[j])
        stuff = stuff + val
        stuff_orig = stuff_orig + val_orig
        j += 1
    invnelec1dm3.append(stuff/np.pi)
    invnelec1dm3_orig.append(stuff_orig/np.pi)
    i += 1

invnelec1dm3 = np.array(invnelec1dm3)
invnelec1dm3_orig = np.array(invnelec1dm3_orig)




# MATRIX TRANSFORMATION SOLUTION
#vector X = matrix A * vector N ----> Find matrix A
print('\n')
A = np.zeros((ny, ny))
deltar = xa[1] - xa[0]
X = np.flip(newbigx, axis=0)
X_orig = np.flip(newbigx_orig, axis=0) 
r = np.flip(xa, axis=0)
i = 0
while i < ny:
    print('\r', i, '/', ny, end='   ')
    rfactor = np.sqrt(r[i]*deltar)
    j = 0
    while j <= i:
        if j == i:
            numfactor = 1
        else:
            numfactor = np.sqrt(2*i + 1 - 2*j) - np.sqrt(2*i - 1 - 2*j)
        A[i,j] = 2*rfactor*numfactor
        j += 1
    i += 1

matA = np.matrix(A)
matX = matrix(X)
matX_orig = matrix(X_orig)
N = matA.I*matX.T
N = np.flip(np.asarray(N), axis=0)
N = N.flatten()
N_orig = matA.I*matX_orig.T
N_orig = np.flip(np.asarray(N_orig), axis=0)
N_orig = N_orig.flatten()


# FIND RESIDUALS
resid_abel = np.abs(nelec2dm3[nc,:] / invnelec1dm3)
resid_abel_orig = np.abs(nelec2dm3_orig[nc,:] / invnelec1dm3_orig)

resid_mat = np.abs(nelec2dm3[nc,:] / N)
resid_mat_orig = np.abs(nelec2dm3_orig[nc,:] / N_orig)

# PLOT
#%%

plt.figure(figsize=(15,9))
plt.subplot(131)
#expected results
plt.semilogx(nelec2dm3_orig[nc,:], yarr1dkm-rpkm, linewidth=2, color='k', alpha=.8, label='expected results (no blob)')
plt.semilogx(nelec2dm3[nc,:], yarr1dkm-rpkm, linewidth=2, color='k', alpha=.8, label='expected results (with blob)')
# Abel transform 
plt.semilogx(invnelec1dm3, yarr1dkm-rpkm, linewidth=5, color='r', alpha=.5, label='Abel transform local density (with blob)')
plt.semilogx(invnelec1dm3_orig, yarr1dkm-rpkm, linewidth=4, color='b', alpha=.5, label='Abel transform local density (no blob)')
# Matrix solution
plt.semilogx(N, yarr1dkm-rpkm, linewidth=3, color='g', alpha=.8, label='Matrix inverse local density (with blob)')
plt.semilogx(N_orig, yarr1dkm-rpkm, linewidth=3, color='orange', alpha=.5, label='Matrix inverse local density (no blob)')
'''
# peaks
plt.hlines(yarr1dkm[invnelec1dm3==np.max(invnelec1dm3)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='r', linestyles='dashed', linewidth=4, label='peak')
plt.hlines(yarr1dkm[invnelec1dm3_orig==np.max(invnelec1dm3_orig)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='b', linestyles='dashed', linewidth=3, label='peak')
plt.hlines(yarr1dkm[np.where(N==np.max(N))[0]]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='g', linestyles='dashed', linewidth=2, label='peak')
plt.hlines(yarr1dkm[np.where(N_orig==np.max(N_orig))[0]]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='orange', linestyles='dashed', label='peak')
'''
plt.xlim(1e-10,1e12)
plt.ylim(0,1000)
plt.legend()
plt.title('Inverse Method Comparison')

plt.subplot(132) # residuals
plt.semilogx(resid_abel, yarr1dkm-rpkm, color='b', label='abel (blob)')

plt.semilogx(resid_mat, yarr1dkm-rpkm, color='g', label='matrix (blob)')
plt.xlim(.5,100)
plt.ylim(0,1000)
plt.legend()
plt.title('Quotient Redisual (blob)')

plt.subplot(133)
plt.semilogx(resid_abel_orig, yarr1dkm-rpkm, color='b', linestyle='dashed', linewidth=3, alpha=.3, label='abel (no blob)')
plt.semilogx(resid_mat_orig, yarr1dkm-rpkm, color='g', linestyle='dashed', linewidth=3, alpha=.3, label='matrix (no blob)')
plt.xlim(.5, 100)
plt.ylim(0,1000)
plt.legend(loc='upper right')
plt.title('Quotient Redisual (no blob)')

#plt.savefig('/Users/saunders/Documents/planet_research/blob_discrete/local_density-scaled_up.png', bbox_inches='tight', dpi=200)
plt.show()

'''
nelec2dm3_b = np.zeros((nx, ny))
rarr2dm = rarr2dkm * 1e3

for i, row in enumerate(rarr2dm):
    interp = np.interp(row, yarr1dm, invnelec1dm3)
    nelec2dm3_b[i,:] = interp

dl = xarr1dm[1] - xarr1dm[0]
tec1dm2_b = np.sum(nelec2dm3_b, axis=0) * dl

plt.figure(figsize=(12,7))
plt.semilogx(tec1dm2_orig, yarr1dkm, 'b', linewidth=7, alpha=.5, label='original')
plt.semilogx(tec1dm2, yarr1dkm, 'r', linewidth=5, alpha=.5, label='blob')
plt.semilogx(tec1dm2_b, yarr1dkm, 'k', linestyle='dashed', label='backwards')
plt.legend()
plt.savefig('/Users/saunders/Documents/planet_research/blob_real/reverse_complicated.png', bbox_inches='tight')
plt.show()
'''