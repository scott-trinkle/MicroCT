function cartop2,fxy, NR = nr, NT = nt, CENTER = center

dims = size(fxy, /dimensions)
if(n_elements(dims) ne 2) then message,"fxy not a 2D array"

if(~keyword_set(nr)) then nr = dims[0]
if(~keyword_set(nt)) then nt = dims[0]

if(~keyword_set(center)) then begin
  x0 = dims[0]/2. - .5
  y0 = dims[1]/2. - .5
endif else begin
  x0 = center[0]
  y0 = center[1]
endelse

r = (findgen(nr)+1)/float(nr)*x0
t = 2.*!pi*findgen(nt)/float(nt-1)


xpolar=x0+r#cos(t)
ypolar=y0+r#sin(t)

;
;
;
; - - since xpolar, ypolar not on a regular xy grid, can't use /grid keyword
; - - in interpolate.
;
frt=interpolate(fxy,xpolar,ypolar,cubic=-0.5)
;
return, frt

end