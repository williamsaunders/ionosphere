### onionpeelviatec.py from onionpeelviatec.pro ###

# imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import ticker, cm
import matplotlib.patches as mpatches
import copy


plt.rcParams['font.size']=14
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

xarr1dkm = (2*np.arange(nx)/(nx-1) - 1) * rpkm 
yarr1dkm = (np.arange(ny)/(ny-1) + 1) * rpkm

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

# get rid of 0 density values so we can take the log
x0, y0 = np.where(nelec2dcm3==0)
nelec2dcm3_0 = copy.copy(nelec2dcm3) #make a copy that we can mess with
nelec2dcm3_0[x0,y0] = 1e-300

# make contour plot of electorn density
#plt.contourf(yarr1dkm, xarr1dkm, nelec2dcm3, np.sort(histedges_equalN(np.concatenate(nelec2dcm3),200)))
#plt.colorbar()
#plt.show()


# make color contour plot for retention
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log(nelec2dcm3_0)
CS = plt.contourf(X, Y, Z, cmap=cm.PuBu_r)
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel(r'log electron density [cm$^{-3}$]')

# formatting and labeling
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
#plt.legend(loc='lower left')
#plt.gca().add_artist(legend1)
#plt.title(r'Disk Evolution in $\alpha$-$\dot{M_{pe}}$ Space')
plt.tight_layout()
plt.show()
plt.savefig('/Users/saunders/Documents/planet_research/contour_sym.png', bbox_inches='tight')



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

invnelec1dm3 = np.array(invnelec1dm3)
plt.semilogx(invnelec1dm3, yarr1dkm-rpkm, linewidth=3)
plt.semilogx(nelec2dm3[nc,:], yarr1dkm-rpkm)
plt.xlim(1e-2,1e12)
plt.ylim(0,400)
plt.show()
plt.savefig('/Users/saunders/Documents/planet_research/coldensity_sym.png', bbox_inches='tight')
