; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

npts = 21
;col = (/"black","red","orange","blue","paleturquoise","green","olivedrab"/)
;dum1 = new(6, graphic)
;dum2 = dum1

season=(/"DJF","JJA"/)
;period = (/"-60","60","-30","30","Inf"/)
;period = (/"_inf","50","20"/)
period = (/"Inf"/)

lats = -ispan(-80, 80, 5)
lons = ispan(60, 360, 90)

years = ispan(1980,1985,1)
nyr = dimsizes(years)

do ssn = 0,0
 do ip    = 0, dimsizes(period)-1
  p = period(ip)
   do lat  = 0,dimsizes(lats)-1
   do lon  = 0,dimsizes(lons)-1
    do k =  1, 6
     fout = "../output/matlab/yearly/raydens_"+ssn+years(0)+"_"+years(nyr-1)+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k
     ; ***  Plotting  *******************************************************************
       wks  = gsn_open_wks("png",fout)

       res                    = True
       res@gsnMaximize        = True         ; make ps, pdf, eps, .. large
       res@gsnDraw            = False        ; don't draw plot yet
       res@gsnFrame           = False        ; don't advance frame yet

       ;pres@vpHeightF= 0.4                    ; change aspect ratio of plot
       ;res@vpWidthF = 5.


       res@mpMinLatF          = 0    ; range to zoom in on
       res@mpMaxLatF          = 90
       res@mpMinLonF          = 0
       res@mpMaxLonF          = 360

       res@tiMainString       = "Ray Path"
       res@tiMainFontHeightF  = 0.013
       ;res@gsnCenterString    = "("+lats(lat)+"N; "+lons(lon)+"E)      period = "+p+"days     root = "+root
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

     do yr = years(0),years(nyr-1)
      do root =  1,3

      fin = "../output/matlab/yearly/ray_"+season(ssn)+yr+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k+"_root"+root
      print(fin)
      if (fileexists(fin)) then
       print("Raypath from "+fin)
       ray = asciiread (fin,(/npts,11/),"float")
       raylon = ray(:,5)
       raylat = ray(:,6)
       print(raylon)
       print(raylat)

       opt          = True
       opt@binx_min= 0
       opt@binx_max= 360
       opt@biny_min= -90
       opt@biny_max=  90
       pdf2 = pdfxy(raylon, raylat, 360, 180, opt)

       print (pdf2@binx_center)
       print (pdf2@biny_center)


       print (dimsizes(pdf2))
        do i = 0,359
        do j=0,179
         if (pdf2(j,i).gt.0) then
          print (i+"   "+j+"   "+pdf2(j,i))
          print (pdf2@binx_center(i)+"  "+pdf2@biny_center(j))
         end if
        end do
       end do

       print(raylon)
       print(raylat)

       ;plot contour plots!!!!
       ; think of a better value for pdf2 7 -> 1
       ; take average of 5 years + all roots




       status_exit(1)

       ; gcdist = gc_latlon(lat1,lon1, lat2,lon2, npts,2)
       ; print (gcdist@gclat+"  "+gcdist@gclon )

       pres@gsLineColor  = col(k-1)         ; color of lines
       pres@xyExplicitLegendLabels = (/"line3","line4"/)
       ; dum1 = gsn_add_polyline(wks,plot, gcdist@gclon ,gcdist@gclat ,pres)
       dum1(k-1) = gsn_add_polyline(wks, plot, raylon, raylat, pres)

       ; Add markers to each of the 'npts' on the great circle path.
       mkres@gsMarkerColor  = col(k-1)
       ; dum2 = gsn_add_polymarker(wks,plot,gcdist@gclon ,gcdist@gclat ,mkres)
       dum2(k-1) = gsn_add_polymarker(wks, plot, raylon, raylat, mkres)

       lres@gsMarkerColor  = col(k-1)
       ;gsn_polymarker_ndc(wks,xleg(k-1),yleg(k-1),lres)
       ;gsn_text_ndc      (wks,"k"+k,xtxt(k-1),ytxt(k-1),txres)
       draw(plot)
       ;frame(wks)
       delete(ray)
       delete(raylon)
       delete(raylat)

       end if  ; if filein exist
      end do  ; root
     end do  ;year

;---Drawing the plot will draw the attached polylines and markers.
  draw(plot)
  frame(wks)

; ***  End plotting  ***************************************************************

    end do  ; k
   end do  ;lon
   end do  ;lat
 end do  ;ip
end do   ; ssn

end
