begin
;************************************************
; read in netCDF file
;************************************************
 var = "v"
  a = addfile("../output/test/test."+var+".nc","r")
;************************************************
  v = a->v
;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks("png","../output/test/test."+var)                  ; send graphics to PNG file

  res = True
  res@mpMinLatF            = -90      ; range to zoom in on
  res@mpMaxLatF            =  90.
  res@mpMinLonF            =  0.
  res@mpMaxLonF            = 360.

  ;res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  ;res@cnMinLevelValF  = -500.			; set the minimum contour level
  ;res@cnMaxLevelValF  = 500.			; set the maximum contour level
  ;res@cnLevelSpacingF = 10.			; set the interval between contours

  ; print(qx)



  plot = gsn_csm_contour_map(wks,v,res)         ; create a default plot
end
