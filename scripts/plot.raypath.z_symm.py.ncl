; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

npts = 31
col = (/"black","red","orange","blue","paleturquoise","green","olivedrab"/)
dum1 = new(6, graphic)
dum2 = dum1

lon = (/"45","180","120"/)
lat = (/"15","20","50"/)
  do loc  = 1,1
   do root = 0,1
    ;period = (/"-60","60","-30","30","Inf"/)
    period = (/"_inf","-14","14"/)
    ;period = (/"-14"/)
    do ip    = 0, dimsizes(period)-1
     p = period(ip)
     fout = "../output/zon_symm/ray_zs_vln_loc"+lat(loc)+"N_"+lon(loc)+"E"+"_period"+p+"_root"+root

; ***  Plotting  *******************************************************************
  wks  = gsn_open_wks("png",fout)

  res                    = True
  res@gsnMaximize        = True         ; make ps, pdf, eps, .. large
  res@gsnDraw            = False        ; don't draw plot yet
  res@gsnFrame           = False        ; don't advance frame yet

  res@mpMinLatF          = -90   ; range to zoom in on
  res@mpMaxLatF          = 90
  res@mpMinLonF          = 0
  res@mpMaxLonF          = 360

  res@tiMainString       = "Ray Path"
  res@tiMainFontHeightF  = 0.013
  res@gsnCenterString    = "("+lat(loc)+"N; "+lon(loc)+"E)      period = "+p+"days     root = "+root
  res@gsnCenterStringFontHeightF = 0.008


  res@mpGreatCircleLinesOn = False

  plot = gsn_csm_map_ce(wks,res)

  ; resources for polylines
  pres = True

  ; resources for markers
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.02


; legend

  lres               = True
  lres@gsMarkerIndex = 17     ; Filled circle
  lres@gsMarkerSizeF = 0.03

  txres               = True
  txres@txFontHeightF = 0.013

  xleg = (/0.15,0.15,0.35,0.35,0.56,0.56/)   ; Location of
  ;xtxt = (/0.22,0.225,0.42,0.44,0.65,0.66/)  ; legend markers
  xtxt = xleg + 0.07
  yleg = (/0.16,0.10,0.16,0.10,0.16,0.10/)   ; and text
  ytxt = yleg
; *** Map plotted ******************************************************************

     do k    = 1,6
     print("k="+k)

      fin = "../output/zon_symm/raypath_zs_vln_loc"+lat(loc)+"N_"+lon(loc)+"E"+"_period"+p+"_k"+k+"_root"+root
      if (fileexists(fin)) then
      print("Raypath from "+fin)
       ray = asciiread (fin,(/npts,29/),"float")
       raylon = ray(:,5)
       raylat = ray(:,6)

       print(raylon)
       print(raylat)


       ; gcdist = gc_latlon(lat1,lon1, lat2,lon2, npts,2)
       ; print (gcdist@gclat+"  "+gcdist@gclon )

       pres@gsLineColor  = col(k-1)         ; color of lines
       pres@xyExplicitLegendLabels = (/"line3","line4"/)
       ; dum1 = gsn_add_polyline(wks,plot, gcdist@gclon ,gcdist@gclat ,pres)
       dum1(k-1) = gsn_add_polyline(wks, plot, raylon, raylat, pres)

       ; Add markers to each of the 'npts' on the great circle path.
       mkres@gsMarkerColor  = "black"
       ; dum2 = gsn_add_polymarker(wks,plot,gcdist@gclon ,gcdist@gclat ,mkres)
       dum2(k-1) = gsn_add_polymarker(wks, plot, raylon, raylat, mkres)

       lres@gsMarkerColor  = col(k-1)
       gsn_polymarker_ndc(wks,xleg(k-1),yleg(k-1),lres)
       gsn_text_ndc      (wks,"k"+k,xtxt(k-1),ytxt(k-1),txres)
       draw(plot)
       ;frame(wks)
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
