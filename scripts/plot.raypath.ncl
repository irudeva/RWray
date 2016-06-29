; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   

begin

npts = 30
col = (/"red","blue","orange","green","black"/)

  do loc  = 2,2
   do root = 1,1
    ;period = (/"-60","60","-30","30","Inf"/)
    period = (/"-60"/)
    do ip    = 0, dimsizes(period)-1
     p = period(ip)
     fout = "../output/plot/ray_location"+loc+"_period"+p+"_root"+root

; ***  Plotting  *******************************************************************
  wks  = gsn_open_wks("png",fout)

  res                    = True
  res@gsnMaximize        = True         ; make ps, pdf, eps, .. large
  res@gsnDraw            = False        ; don't draw plot yet
  res@gsnFrame           = False        ; don't advance frame yet

  res@mpMinLatF          = 20   ; range to zoom in on
  res@mpMaxLatF          = 90
  res@mpMinLonF          = -180
  res@mpMaxLonF          = 180

  res@tiMainString       = "Ray Path: period = "+p+"days location = "+loc+" root = "+root
  res@tiMainFontHeightF  = 0.015
  ;res@gsnCenterString    = "res@mpGreatCircleLinesOn = True"
  ;res@gsnCenterStringFontHeightF = 0.011


  res@mpGreatCircleLinesOn = False

  plot = gsn_csm_map_ce(wks,res)  

  ; resources for polylines
  pres = True  
  ; resources for merkers
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.035
; *** Map plotted ******************************************************************

     do k    = 1,1
    
      fin = "../output/raypath_k"+k+"_period"+p+"_location"+loc+"_root"+root
      print("Raypath from "+fin)
      if (fileexists(fin)) then
       ray = asciiread (fin,(/30,3/),"float")
       raylon = ray(:,1)
       raylat = ray(:,2)

       ; gcdist = gc_latlon(lat1,lon1, lat2,lon2, npts,2)
       ; print (gcdist@gclat+"  "+gcdist@gclon ) 
  
       pres@gsLineColor  = col(k)         ; color of lines
       ; dum1 = gsn_add_polyline(wks,plot, gcdist@gclon ,gcdist@gclat ,pres)
       dum1 = gsn_add_polyline(wks, plot, raylon, raylat, pres)

       ; Add markers to each of the 'npts' on the great circle path.
       mkres@gsMarkerColor  = "black" 
       ; dum2 = gsn_add_polymarker(wks,plot,gcdist@gclon ,gcdist@gclat ,mkres)
       dum2 = gsn_add_polymarker(wks, plot, raylon, raylat, mkres)
 
       delete(ray)
       delete(raylon)
       delete(raylat)

     end if  ; if filein exist
    end do  ; k

;---Drawing the plot will draw the attached polylines and markers.
  draw(plot)
  frame(wks)

; ***  End plotting  ***************************************************************

   end do  ; ip
  end do  ;root
 end do  ;loc

end
