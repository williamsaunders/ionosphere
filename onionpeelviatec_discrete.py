### onionpeelviatec.py from onionpeelviatec.pro ###

# imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import ticker, cm
import matplotlib.patches as mpatches
import copy
from numpy import matrix


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

# convert to m-3
nelec2dm3 = nelec2dcm3 * 1e6
n0m3 = n0cm3 * 1e6

# get rid of 0 density values so we can take the log
x0, y0 = np.where(nelec2dm3==0)
nelec2dm3_0 = copy.copy(nelec2dm3) #make a copy that we can mess with
nelec2dm3_0[x0,y0] = 1e-300

# make contour plot of electorn density
plt.figure(figsize=(12,7))
X = xarr2dkm
Y = yarr2dkm
Z = np.log10(nelec2dm3_0)
Z2 = np.log10(n0m3) + (-1 - dummyx2d - np.exp(-dummyx2d))*np.log10(np.e)
contours = [np.min(Z2)] + np.arange(-100,30,10).tolist()
CS = plt.contourf(X, Y, Z2,  contours, cmap=cm.PuBu_r, vmin=-100, vmax=20)
cbar = plt.colorbar(CS,ticks=contours)
cbar.ax.set_ylabel(r'log electron density [m$^{-3}$]')

# formatting and labeling
plt.xlabel('X position [km]', fontsize=20)
plt.ylabel('Y position [km]', fontsize=20)
plt.gcf().subplots_adjust(bottom=0.15)
plt.grid()
plt.tight_layout()
#plt.savefig('/Users/saunders/Documents/planet_research/contour_sym-whole.png', bbox_inches='tight', dpi=200)
plt.ylim(3400,3800)
#plt.savefig('/Users/saunders/Documents/planet_research/contour_sym-low.png', bbox_inches='tight', dpi=200)
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

# vector X = matrix A * vector N ----> Find matrix A

A = np.zeros((ny, ny))
deltar = xa[1] - xa[0]
X = np.flip(newbigx, axis=0)
r = np.flip(xa, axis=0)
i = 0
while i < ny:
    print(i, '/', ny)
    rfactor = np.sqrt(r[i]*deltar)
    j = 0
    while j <= i:
        if j == i:
            numfactor = 1
        else:
            numfactor = np.sqrt(2*i + 1 - 2*j) - np.sqrt(2*i - 1 - 2*j)
        A[i,j] = rfactor*numfactor
        j += 1
    i += 1

matA = np.matrix(A)
matX = matrix(X)
N = matA.I*matX.T
N = np.asarray(N)
N = np.flip(N, axis=0)



# perform integration

invnelec1dm3 = []
i = 0
while i < len(newrevdxdr):
    stuff = 0.
    print(i, '/', len(newrevdxdr))
    j = i
    while j < len(newrevdxdr) - 1:
        stuff = stuff + .5 * ( np.log(reva2[j]/reva2[i] + np.sqrt((reva2[j]/reva2[i])**2. - 1.)) + \
                np.log(reva2[j+1]/reva2[i] + np.sqrt((reva2[j+1]/reva2[i])**2. - 1.)) ) \
                * (newrevdxdr[j+1] - newrevdxdr[j])
        j += 1
    invnelec1dm3.append(stuff/np.pi)
    i += 1

invnelec1dm3 = np.array(invnelec1dm3)



plt.figure(figsize=(12,7))
plt.semilogx(invnelec1dm3, yarr1dkm-rpkm, linewidth=3, color='k')
plt.semilogx(nelec2dm3[nc,:], yarr1dkm-rpkm, color='b')
plt.semilogx(N, yarr1dkm-rpkm, linewidth=5, alpha=.5, color='g')
plt.xlim(1e-2,1e12)
plt.ylim(0,400)
plt.savefig('/Users/saunders/Documents/planet_research/matrix_density.png', bbox_inches='tight')
plt.show()

