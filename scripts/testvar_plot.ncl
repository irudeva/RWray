begin
;************************************************
; read in netCDF file
;************************************************
 var = "umx"
  a = addfile("../output/test/test."+var+".nc","r")
;************************************************
  umx = a->umx
;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks("png","../output/test/test."+var)                  ; send graphics to PNG file

  res = True
  res@mpMinLatF            = -90      ; range to zoom in on
  res@mpMaxLatF            =  90.
  res@mpMinLonF            =  0.
  res@mpMaxLonF            = 360.

  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = -20.			; set the minimum contour level
  res@cnMaxLevelValF  = 20.			; set the maximum contour level
  res@cnLevelSpacingF = 2.			; set the interval between contours



  plot = gsn_csm_contour_map(wks,umx,res)         ; create a default plot
end
