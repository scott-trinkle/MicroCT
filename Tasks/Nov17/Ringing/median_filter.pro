function median_filter, array, width 
on_error, 1

if(n_elements(width) eq 1) then return, median(array, width)
if(size(width, /n_dimensions) gt 1) then message, 'Width dimesion must be 1'

return, 1

end
