### onionpeelviatec.py from onionpeelviatec.pro ###

# imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import ticker, cm
import matplotlib.patches as mpatches
import copy

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
theta2ddeg = theta2rad * np.pi/180.

# solar zenith angle not used at the moment 

# define spherical ionosphere                                                                                                         
dummyx2d = (rarr2dkm - r0km) / dsh0km
nelec2dcm3 = n0cm3 * np.exp(1 - dummyx2d - np.exp(-dummyx2d))

# add blob of plasma modeled as a gaussian distribution
sigma = .03*nx
#mu = (np.random.choice(range(nx)), np.random.choice(ny))
#mu = (int(nx/2.),ny-int(ny/10.))
mu = (200,200)
xg = (1/np.sqrt(2*np.pi)) * np.exp(-.5*((xarr1dkm-xarr1dkm[mu[0]])/sigma)**2)
yg = (1/np.sqrt(2*np.pi)) * np.exp(-.5*((yarr1dkm-yarr1dkm[mu[1]])/sigma)**2)
Xg, Yg= np.meshgrid(yg,xg)
peak = nelec2dcm3[mu[0],mu[1]]*1e5
peak2 = np.max(nelec2dcm3)*1e-25
nelec2dcm3_g = n0cm3*Xg*Yg*peak2
nelec2dcm3_orig = copy.copy(nelec2dcm3)
nelec2dcm3 = nelec2dcm3 + nelec2dcm3_g
#nelec2dcm3 = nelec2dcm3_g

# convert to m-3 
nelec2dm3 = nelec2dcm3 * 1e6
nelec2dm3_g = nelec2dcm3_g * 1e6
nelec2dm3_orig = nelec2dcm3_orig*1e6
n0m3 = n0cm3 * 1e6

'''
# get rid of 0 density values so we can take the log
x0, y0 = np.where(nelec2dm3_g==0)
nelec2dm3_g0 = copy.copy(nelec2dm3_g)
nelec2dm3_g[x0,y0] = 1e-300
x1, y1 = np.where(nelec2dm3==0)
nelec2dm3_0 = copy.copy(nelec2dm3) 
nelec2dm3_0[x1,y1] = 1e-300
'''

# plot just the blob
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log10(nelec2dcm3_g)
x0, y0 = np.where(Z==-np.inf)
Z2 = copy.copy(Z)
Z2[x0,y0] = -70000
#contours = [np.min(Z2)] + np.arange(-100,30,10).tolist()
contours = [-70000] + np.arange(-100,30,10).tolist()
CS = plt.contourf(X, Y, Z2,  contours, cmap=cm.PuBu_r, vmin=-100, vmax=20)
cbar = plt.colorbar(CS,ticks=contours)
cbar.ax.set_ylabel(r'log electron density [cm$^{-3}$]')
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
plt.tight_layout()
#plt.savefig('/Users/saunders/Documents/planet_research/blob.png', bbox_inches='tight')
plt.show()

# make contour plot of electorn density
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log10(nelec2dm3)
x0, y0 = np.where(Z==-np.inf)
Z2 = copy.copy(Z)
Z2[x0,y0] = -70000
#Z2 = np.log10(n0m3) + (-1 - dummyx2d - np.exp(-dummyx2d))*np.log10(np.e)
#contours = [np.min(Z2)] + np.arange(-100,30,10).tolist()
CS = plt.contourf(X, Y, Z2,  contours, cmap=cm.PuBu_r, vmin=-100, vmax=20)
cbar = plt.colorbar(CS,ticks=contours)
cbar.ax.set_ylabel(r'log electron density [m$^{-3}$]')

# formatting and labeling
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
plt.tight_layout()
plt.savefig('/Users/saunders/Documents/planet_research/blob_test/med_blob.png', bbox_inches='tight', dpi=200)
#plt.ylim(3400,3800)
#plt.savefig('/Users/saunders/Documents/planet_research/contour_blob-low.png', bbox_inches='tight', dpi=200)
plt.show()


# unit conversions
xarr1dm = xarr1dkm * 1e3
yarr1dm = yarr1dkm * 1e3

dl = xarr1dm[1] - xarr1dm[0]
tec1dm2 = np.sum(nelec2dm3, axis=0) * dl

# change variable names 
newbigx = np.array(tec1dm2)
xa = np.array(yarr1dm)

# perform derivative: dOmega/dX
newdxdr = (np.roll(newbigx,-1) - np.roll(newbigx,1)) / (np.roll(xa,-1) - np.roll(xa,1))
newdxdr[0] = 2.*newdxdr[2] - newdxdr[1]
newdxdr[-1] = 2.*newdxdr[-3] - newdxdr[-2]

# some math checks 
pp = np.argsort(xa) # make sure order is correct becuase sometimes radio occultations \ 
                    # come out of order becuase of refraction or direction of transit 
reva2 = xa[pp] # ordered impact parameter, m (closest approach distance)
thisrkm = reva2 / 1e3 # ordered impact parameter, km
newrevdxdr = newdxdr[pp] # ordered spatial derivative of bigx, m-3

# set first element of newrevdxdr to 0 to avoid problems with numerical integration
newrevdxdr[-1] = 0.

# electrondensity is the local election density m-3, from corrected data
# perform integration

invnelec1dm3 = []
i = 0
while i < len(newrevdxdr):
    stuff = 0.
#    print(i, '/', len(newrevdxdr))
    j = i
    while j < len(newrevdxdr) - 1:
        stuff = stuff + .5 * ( np.log(reva2[j]/reva2[i] + np.sqrt((reva2[j]/reva2[i])**2. - 1.)) + \
                np.log(reva2[j+1]/reva2[i] + np.sqrt((reva2[j+1]/reva2[i])**2. - 1.)) ) \
                * (newrevdxdr[j+1] - newrevdxdr[j])
        j += 1
    invnelec1dm3.append(stuff/np.pi)
    i += 1

# DO EVERYTHING AGAIN FOR THE ORIGINAL  (NO BLOB)

# unit conversions 
nelec2dm3 = nelec2dcm3_orig * 1e6
xarr1dm = xarr1dkm * 1e3
yarr1dm = yarr1dkm * 1e3

dl = xarr1dm[1] - xarr1dm[0]
tec1dm2 = np.sum(nelec2dm3, axis=0) * dl

# change variable names 
newbigx = np.array(tec1dm2)
xa = np.array(yarr1dm)

# perform derivative: dOmega/dX
newdxdr = (np.roll(newbigx,-1) - np.roll(newbigx,1)) / (np.roll(xa,-1) - np.roll(xa,1))
newdxdr[0] = 2.*newdxdr[2] - newdxdr[1]
newdxdr[-1] = 2.*newdxdr[-3] - newdxdr[-2]

# some math checks 
pp = np.argsort(xa) # make sure order is correct becuase sometimes radio occultations \ 
                    # come out of order becuase of refraction or direction of transit 
reva2 = xa[pp] # ordered impact parameter, m (closest approach distance)
thisrkm = reva2 / 1e3 # ordered impact parameter, km
newrevdxdr = newdxdr[pp] # ordered spatial derivative of bigx, m-3

# set first element of newrevdxdr to 0 to avoid problems with numerical integration
newrevdxdr[-1] = 0.

# electrondensity is the local election density m-3, from corrected data

# perform integration
invnelec1dm3_orig = []
i = 0
while i < len(newrevdxdr):
    stuff = 0.
    j = i
    while j < len(newrevdxdr) - 1:
        stuff = stuff + .5 * ( np.log(reva2[j]/reva2[i] + np.sqrt((reva2[j]/reva2[i])**2. - 1.)) + \
                np.log(reva2[j+1]/reva2[i] + np.sqrt((reva2[j+1]/reva2[i])**2. - 1.)) ) \
                * (newrevdxdr[j+1] - newrevdxdr[j])
        j += 1
    invnelec1dm3_orig.append(stuff/np.pi)
    i += 1

invnelec1dm3 = np.array(invnelec1dm3)

plt.figure(figsize=(10,7))
plt.semilogx(invnelec1dm3, yarr1dkm-rpkm, linewidth=5, color='r', label='blob integrated electron column density')
plt.hlines(yarr1dkm[invnelec1dm3==np.max(invnelec1dm3)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='r', linestyles='\
dashed', linewidth=3, label='peak electron density')
plt.semilogx(invnelec1dm3_orig, yarr1dkm-rpkm, linewidth=3, color='b', label='original integrated electron column density')
plt.hlines(yarr1dkm[invnelec1dm3_orig==np.max(invnelec1dm3_orig)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='b', linestyles='dashed', label='original peak electron density')
plt.semilogx(nelec2dm3_orig[nc,:], yarr1dkm-rpkm, linewidth=1, color='g', label='predicted electron column density')                              
plt.xlim(1e-14,1e12)
plt.ylim(0,1000)
plt.legend()
plt.savefig('/Users/saunders/Documents/planet_research/blob_test/med-blob.png', bbox_inches='tight')
plt.show()

