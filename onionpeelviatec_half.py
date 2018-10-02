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
nelec2dcm3_orig = copy.copy(nelec2dcm3)

# make right half higher density 
nelec2dcm3[int(nx/2):,:] = nelec2dcm3_orig[int(nx/2):,:]*1e5

# get rid of 0 density values so we can take the log
x0, y0 = np.where(nelec2dcm3==0)
nelec2dcm3_0 = copy.copy(nelec2dcm3) #make a copy that we can mess with
nelec2dcm3_0[x0,y0] = 1e-300

# make contour plot of electorn density
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log(nelec2dcm3_0)
CS = plt.contourf(X, Y, Z, [-400,-350,-300,-250,-200,-150,-100,-50,0,50,100], cmap=cm.PuBu_r, vmin=-400, vmax=100)
cbar = plt.colorbar(CS,ticks=[-400,-350,-300,-250,-200,-150,-100,-50,0,50,100])
cbar.ax.set_ylabel(r'log electron density [cm$^{-3}$]')

# formatting and labeling
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
#plt.legend(loc='lower left')
#plt.gca().add_artist(legend1)
plt.tight_layout()
plt.savefig('/Users/saunders/Documents/planet_research/contour_half2.png', bbox_inches='tight')
plt.show()


# unit conversions 
nelec2dm3 = nelec2dcm3 * 1e6
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
plt.semilogx(invnelec1dm3, yarr1dkm-rpkm, linewidth=3, color='r', label='blob integrated electron column density')
plt.hlines(yarr1dkm[invnelec1dm3==np.max(invnelec1dm3)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='r', linestyles='dashed', linewidth=3, label='peak electron density')
plt.semilogx(invnelec1dm3_orig, yarr1dkm-rpkm, linewidth=1, color='b', label='original integrated electron column density')
plt.hlines(yarr1dkm[invnelec1dm3_orig==np.max(invnelec1dm3_orig)]-rpkm, np.min(invnelec1dm3), np.max(invnelec1dm3),color='b', linestyles='dashed', label='original peak electron density')
#plt.semilogx(nelec2dm3[nc,:], yarr1dkm-rpkm, color='b', label='predicted electron column density')
plt.xlim(1e-2,1e16)
plt.ylim(0,400)
plt.legend()
plt.savefig('/Users/saunders/Documents/planet_research/elec_half2.png', bbox_inches='tight')
plt.show()
