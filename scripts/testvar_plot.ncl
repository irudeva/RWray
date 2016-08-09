begin
;************************************************
; read in netCDF file
;************************************************
 var = "umy1"
  a = addfile("../output/test/test."+var+".nc","r")
;************************************************
  umy1 = a->umy1                               ; read July zonal winds
;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks("png","../output/test/test."+var)                  ; send graphics to PNG file

  res = True
  res@mpMinLatF            = -90      ; range to zoom in on
  res@mpMaxLatF            =  90.
  res@mpMinLonF            =  0.
  res@mpMaxLonF            = 360.


  plot = gsn_csm_contour_map(wks,umy1,res)         ; create a default plot
end
