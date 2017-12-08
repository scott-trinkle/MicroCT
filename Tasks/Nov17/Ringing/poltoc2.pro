function poltoc2, frt, NX=nx, NY=ny, CENTER=center

dims = size(frt, /dimensions)
if(n_elements(dims) ne 2) then message,"frt not a 2D array"

if(~keyword_set(nx)) then nx = dims[0]
if(~keyword_set(ny)) then ny = dims[0]

if(~keyword_set(center)) then begin
  x0 = nx/2. - .5
  y0 = nx/2. - .5
endif else begin
  x0 = center[0]
  y0 = center[1]
endelse

nr = dims[0]
nt = dims[1]

r = (findgen(nr)+1)/float(nr)*nx
t = 2.*!pi*findgen(nt)/float(nt-1)

xpolar=r#cos(t)
ypolar=r#sin(t)
triangulate,xpolar,ypolar,triangles
fxy=trigrid(xpolar,ypolar,frt,triangles,nx=nx,ny=ny)

return, fxy
end