pro onionpeelviatec

; Paul Withers, 2018.09.10
; Astronomy Department, Boston University

; Explore effects on horizontal gradients on results of Abel transform inversion

; Likely bugs:
; Current formulation of Abel transform (onion-peeling algorithm)
; fails when dNe/dz is large and positive (bottomside of ionosphere)
; See discussion near Equation C13 of Fjeldbo, Kliore, and Eshleman (1971)
; Astronomical Journal, 76, 123-140

; Useful papers (for the old, old ray-tracing version):
; Fjeldbo, Kliore, and Eshleman (1971)
; Astronomical Journal, 76, 123-140
;
; Hinson et al. (1999) JGR, 104, p26997-27012
; 
; Reilly and Strobel (1988)
; Radio Science, 23, 247-256
;
; Kelso (1964) Radio Ray Propagation in the Ionosphere, McGraw-Hill, p260 or p620
;
; Born and Wolf (1965) Optics, 3rd ed., QC355 F.651 in Sci/Eng Library, p733



rpkm = double(3400.) ; Rmars in km
r0km = 120. + rpkm
; altitude of the electron density peak at the subsolar point, km
; converted into radial distance, km
; Hantsch and Bauer (1990) PSS 38, 539-542, figure 3
n0cm3 = 2e5
; peak electron density at the subsolar point, cm-3
; Hantsch and Bauer (1990) PSS 38, 539-542, figure 1
dsh0km = double(10.) ; Density scale height, km
; Hantsch and Bauer (1990) PSS 38, 539-542, figure 3

; Use a Chapman layer ionosphere for simplicity

nx = 2001d
nc = round((nx-1d)/2) ; the central value
ny = 1001d

xarr1dkm = (2d*dindgen(nx)/(nx-1d) - 1d) * rpkm ; -2 to +2 -3400 to 3400 
yarr1dkm = (dindgen(ny)/(ny-1d) + 1d) * rpkm ; 1 to 2       3400 to 6800

xarr2dkm = dblarr(nx,ny)
yarr2dkm = xarr2dkm

i=0
while i lt nx do begin

xarr2dkm[i,*] = xarr1dkm[i]
yarr2dkm[i,*] = yarr1dkm[*]
i++
endwhile

rarr2dkm = sqrt(xarr2dkm^2 + yarr2dkm^2)
theta2drad = atan(yarr2dkm, xarr2dkm)
theta2ddeg = theta2drad * !pi/180d

; not used at the moment
szaoccdeg = 90d ; solar zenith angle at occultation point (x=0, y=rp)
sza2ddeg = abs(90d - theta2ddeg + szaoccdeg)
; this expression works for some range of szaoccdeg around 90, but is not general, fix if you can

if min(sza2ddeg) lt 0. or max(sza2ddeg) gt 180. then stop
; if stop here, then fix the problem - SZA must be between 0 and 180 degrees

dummyx2d = (rarr2dkm - r0km) / dsh0km
nelec2dcm3 = n0cm3 * exp(1d - dummyx2d - exp(-1d * dummyx2d))

nelec2dm3 = nelec2dcm3 * 1d6
xarr1dm = xarr1dkm * 1d3
yarr1dm = yarr1dkm * 1d3

dl = xarr1dm[1] - xarr1dm[0]
tec1dm2 = total(nelec2dm3,1) * dl

; now introducing a bunch of new variable names in order to incorporate some existing code more easily

newbigx = tec1dm2
xa = yarr1dm

newdxdr = (shift(newbigx,-1) - shift(newbigx,1)) / (shift(xa,-1)-shift(xa,1))
newdxdr(0) = 2.*newdxdr(2)-newdxdr(1)
newdxdr(n_elements(newdxdr)-1) = 2.*newdxdr(n_elements(newdxdr)-3) - newdxdr(n_elements(newdxdr)-2)

pp = (sort(xa))
reva2 = xa(pp) ; Ordered impact parameter, m
thisrkm = reva2 / 1e3 ; Ordered impact parameter, km
newrevdxdr = newdxdr(pp) ; Ordered spatial derivative of bigx, m-3

; Set first element of newrevdxdr to zero to avoid problems upon numerical integration
newrevdxdr(n_elements(newrevdxdr)-1) = 0.

; electrondensity is the local electron density, m-3, from corrected data
invnelec1dm3 = newrevdxdr*0.
i=0
while i lt n_elements(invnelec1dm3) do begin
stuff=0.d
j=i
while j lt n_elements(invnelec1dm3)-1 do begin
; Results are very sensitive to details of integration routine
stuff = stuff + 0.5 * ( alog(reva2[j]/reva2[i] + sqrt((reva2[j]/reva2[i])^(2.d)-(1.d)) ) + $
 alog(reva2[j+1]/reva2[i] + sqrt((reva2[j+1]/reva2[i])^(2.d)-(1.d)) )) * (newrevdxdr[j+1]-newrevdxdr[j]) 
j++
endwhile
invnelec1dm3[i] = stuff / !pi
i++
endwhile

plot, invnelec1dm3, yarr1dkm-rpkm, yra=[0,400], thick=3, /xlog
oplot, nelec2dm3[nc,*], yarr1dkm-rpkm, color=255, thick=1
; they agree!


stop
end
