function de_ring, f, w_rd, w_az, NR = nr, NT = nt, $
  THRESHOLD_INITIAL = threshold_inital, THRESHOLD_ARTIFACT = threshold_artifact, THRESH_ZERO = thresh_zero, $
  CENTER = center, STEP = step, VAR = var, $
  R_MED = r_med, RA_MED = ra_med, ART = art
    
;on_error, 2

;Threshold inital image
if(~keyword_set(threshold_initial)) then begin
  threshold_initial    = fltarr(2) 
  threshold_initial[0] = min(f)
  threshold_initial[1] = max(f)
endif
fs = f > threshold_initial[0] < threshold_initial[1]

dims = size(f, /dimensions)


;Polar Transform
frt = cartop2(fs, NR = nr, NT = nt, CENTER = center)

;Radial Median Filter
med_r = fltarr(nr, nt)

if(keyword_set(step)) then begin
  a = fltarr(nr, nt)
  b = fltarr(nr, nt)
  c = fltarr(nr, nt)

  for i=0, nt-1 do a[*,i] =  median(reform(frt[*,i]), (1/3.)*w_rd)
  for i=0, nt-1 do b[*,i] =  median(reform(frt[*,i]), (2/3.)*w_rd)
  for i=0, nt-1 do c[*,i] =  median(reform(frt[*,i]), (3/3.)*w_rd)
  med_r[0:nr/3-1,*]      = a[0:nr/3-1,*]
  med_r[nr/3:2*nr/3-1,*] = b[nr/3:2*nr/3-1,*]
  med_r[2*nr/3:nr-1,*]   = c[2*nr/3:nr-1,*]
endif else if(keyword_set(var)) then begin
  for i =0, nt -1 do med_r[*,i] = median_filter_var(reform(frt[*,i]), [w_rd/3, w_rd])
endif else $
  for i =0, nt -1 do med_r[*,i] = median(reform(frt[*,i]), w_rd)

;Threshold artifact image
art = frt - med_r

if(keyword_set(threshold_artifact)) then begin
  if(keyword_set(thresh_zero)) then begin
    pos = where(art le threshold_artifact[0] or art ge threshold_artifact[1])
    art[pos] = 0.0
  endif else art = art >threshold_artifact[0] <threshold_artifact[1]
endif

;Azimutah Median Filter (Low-pass Filter)
med_ra = fltarr(nr, nt)
if(keyword_set(step)) then begin
  for i = 0, nr/3 - 1 do med_ra[i,*]      = median(reform(art[i,*]), (1/3.)*w_az)
  for i = nr/3, 2*nr/3 - 1 do med_ra[i,*] = median(reform(art[i,*]), (2/3.)*w_az)
  for i = 2*nr/3, nr - 1 do med_ra[i,*]   = median(reform(art[i,*]), (3/3.)*w_az)
endif else if(keyword_set(var)) then begin
  filter = indgen(w_az - w_az/3 + 1) + w_az/3
  nf = n_elements(filter)
  
  for i =0, nr -1 do begin
     j = round(float(i)/(nr-1)*(nf-1))
     large_med_tmp = median(shift([reform(art[i,*]), reform(art[i,*])], nt/2 ), filter[j])
     large_med_tmp = shift(large_med_tmp, -nt/2)
     
     med_ra[i,*] = large_med_tmp[0:nt-1]
  endfor
endif else $
  for i =0, nr -1 do med_ra[i,*] = median(reform(art[i,*]), w_az)

med_ra =  poltoc2(med_ra, nx = dims[0], ny = dims[1], CENTER = center)

r_med  = med_r
ra_med = med_ra

return, f - med_ra
end
