; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   

begin
  npts = 10              ; arbitrary
  lat1 = 20.
  lon1 = -120.
  lat2 = 60.
  lon2 = -64.
  
  wks  = gsn_open_wks("png","../output/plot/raypath")    ; send graphics to PNG file 

  res                    = True
  res@gsnMaximize        = True         ; make ps, pdf, eps, .. large
  res@gsnDraw            = False        ; don't draw plot yet
  res@gsnFrame           = False        ; don't advance frame yet

  space                  =   5          ; arbitrary
  res@mpMinLatF          = lat1-space   ; range to zoom in on
  res@mpMaxLatF          = lat2+space
  res@mpMinLonF          = lon1-space
  res@mpMaxLonF          = lon2+space

  res@tiMainString       = "1st method: Two Points and Great Circle Path"
  res@tiMainFontHeightF  = 0.015
  res@gsnCenterString    = "res@mpGreatCircleLinesOn = True"
  res@gsnCenterStringFontHeightF = 0.011


  res@tiMainString       = "2nd method: Two Points and Great Circle Path"
  res@gsnCenterString    = "gc_latlon used to calculate great circle points"

  res@mpGreatCircleLinesOn = False

  plot2 = gsn_csm_map_ce(wks,res)  

  gcdist = gc_latlon(lat1,lon1, lat2,lon2, npts,2)

  print (gcdist)
  print (gcdist@gclat+"  "+gcdist@gclon )  ; print the lats/lons
  
  pres                  = True  
  pres@gsLineColor      = "red"         ; color of lines
  dum2 = gsn_add_polyline(wks,plot2, gcdist@gclon ,gcdist@gclat ,pres)

;
; Additionally add markers to each of the 'npts' on the 
; great circle path.
;
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.035
  mkres@gsMarkerColor  = "black" 
  dum3 = gsn_add_polymarker(wks,plot2,gcdist@gclon ,gcdist@gclat ,mkres)

;---Drawing the plot will draw the attached polylines and markers.
  draw(plot2)
  frame(wks)

end
