; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   

begin

npts = 10
col = (/"red","blue","orange","green","black","black","black"/)

  do loc  = 3,3
   do root = 1,1
    ;period = (/"-60","60","-30","30","Inf"/)
    period = (/"-14","14"/)
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

  res@tiMainString       = "Ray Path"
  res@tiMainFontHeightF  = 0.013
  res@gsnCenterString    = "location "+loc+"      period = "+p+"days     root = "+root
  res@gsnCenterStringFontHeightF = 0.008


  res@mpGreatCircleLinesOn = False

  plot = gsn_csm_map_ce(wks,res)  

  ; resources for polylines
  pres = True  
  ; resources for merkers
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.02
; *** Map plotted ******************************************************************

     do k    = 1,6
    
      fin = "../output/raypath_k"+k+"_period"+p+"_loc"+loc+"_root"+root
      print("Raypath from "+fin)
      if (fileexists(fin)) then
       ray = asciiread (fin,(/npts,3/),"float")
       raylon = ray(:,2)
       raylat = ray(:,1)

       ; gcdist = gc_latlon(lat1,lon1, lat2,lon2, npts,2)
       ; print (gcdist@gclat+"  "+gcdist@gclon ) 
  
       pres@gsLineColor  = col(k-1)         ; color of lines
       ; dum1 = gsn_add_polyline(wks,plot, gcdist@gclon ,gcdist@gclat ,pres)
       dum1 = gsn_add_polyline(wks, plot, raylon, raylat, pres)

       ; Add markers to each of the 'npts' on the great circle path.
       mkres@gsMarkerColor  = "black" 
       ; dum2 = gsn_add_polymarker(wks,plot,gcdist@gclon ,gcdist@gclat ,mkres)
       dum2 = gsn_add_polymarker(wks, plot, raylon, raylat, mkres)
 
       draw(plot)
       frame(wks)
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
