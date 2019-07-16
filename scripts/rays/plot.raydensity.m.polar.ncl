; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

npts = 21
col = (/"black","red","orange","blue","hotpink","green","olivedrab","cadetblue4","coral4","cyan1"/)
;dum2 = dum1

season=(/"DJF","JJA"/)
;period = (/"-60","60","-30","30","Inf"/)
;period = (/"_inf","50","20"/)
period = (/"Inf"/)

;lats = -ispan(-80, 80, 5)
;lons = ispan(60, 360, 90)
lats = (/72,72/)
lons = (/40,160/)

fyr = 1980
lyr = 2017

; do iy = 1981,2016,100

;years = ispan(iy,iy+9,1)
years = ispan(iy,iy,1)
nyr = dimsizes(years)

longrid = fspan(0.,359.,360)
latgrid = fspan(-90,90,181)

roots = ispan(1,3,1)

;onto grid parameters
ngc0 = 10 ;number of points for interpolation (gcdist)
opt          = True
opt@binx_min = -.5
opt@binx_max = 359.5
opt@biny_min = -90.5
opt@biny_max =  90.5

rayfreq = new((/181,360/),integer)
rayfreq!0 = "lat"
rayfreq!1 = "lon"

rayfreq&lat=latgrid
rayfreq&lon=longrid


rayfreq&lat@units="degrees_north"
rayfreq&lon@units="degrees_east"

;dum = new(dimsizes(roots)*dimsizes(years), graphic)



do ssn = 0,0
 do ip    = 0, dimsizes(period)-1
  p = period(ip)
   do lat  = 0,dimsizes(lats)-1
   do lon  = 0,dimsizes(lons)-1
    do k =  1, 6
      ; fout = "../output/matlab/yearly1/raydensP_"+season(ssn)+years(0)+"_"+years(nyr-1)+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k
      fout = "../output/matlab/clim/rayplot/raydensP_"+season(ssn)+years(0)+"_"+years(nyr-1)+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k
     rayfreq = 0
     dum = new(dimsizes(roots)*dimsizes(years), graphic)
     mark = dum

     ; ***  Plotting  *******************************************************************
       wks  = gsn_open_wks("png",fout)

       res                    = True
       ;res@gsnMaximize       = True         ; make ps, pdf, eps, .. large
       res@gsnDraw            = False        ; don't draw plot yet
       res@gsnFrame           = False        ; don't advance frame yet

       ;pres@vpHeightF= 0.4                    ; change aspect ratio of plot
       ;res@vpWidthF = 5.

       res@gsnPolar   = "NH"

       res@mpMinLatF          = 45    ; range to zoom in on
       res@mpMaxLatF          = 90
       res@mpMinLonF          = 0
       res@mpMaxLonF          = 360
       res@mpCenterLonF       = 180

       res@tiMainString       = "Ray Path"
       res@tiMainFontHeightF  = 0.013
       ;res@gsnCenterString    = "("+lats(lat)+"N; "+lons(lon)+"E)      period = "+p+"days     root = "+root
       ;res@gsnCenterStringFontHeightF = 0.008


       res@mpGreatCircleLinesOn = False

       plot = gsn_csm_map_polar(wks,res)

       ; resources for polylines
       pres = True

       ; resources for markers
       mkres               = True
       mkres@gsMarkerIndex = 17     ; Filled circle
       mkres@gsMarkerSizeF = 0.02


      ;legend

       lres               = True
       lres@gsMarkerIndex = 17     ; Filled circle
       lres@gsMarkerSizeF = 0.03

       txres               = True
       txres@txFontHeightF = 0.013

       xleg = (/0.1,0.1,0.25,0.25,0.40,0.40,.55,.55,.7,.7/)   ; Location of
       ;;xtxt = (/0.22,0.225,0.42,0.44,0.65,0.66/)  ; legend markers
       xtxt = xleg + 0.07
       yleg = (/0.16,0.10,0.16,0.10,0.16,0.10,0.16,0.10,0.16,0.10/)   ; and text
       ytxt = yleg
     ; *** Map plotted ******************************************************************

     idum = -1
     do yr = years(0),years(nyr-1)
      do root =  roots(0),roots(dimsizes(roots)-1)

      fin = "../output/matlab/clim/rays/ray_"+season(ssn)+fyr+"-"+lyr+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k+"_root"+root
      print(fin)
      if (fileexists(fin)) then
       print("Raypath from "+fin)
       ray = asciiread (fin,(/npts,15/),"float")
       raylon = ray(:,5)
       raylat = ray(:,6)
       ;print(raylon)
       ;print(raylat)

       ;pres@gsLineColor  = col(yr-years(0))         ; color of lines"
       ;pres@xyExplicitLegendLabels = (/"line3","line4"/)
       ;if (yr.ge.1990.and.yr.lt.2000) then
       ;   pres@gsLineColor = col(yr-years(0)-10)         ; color of lines"
      ;   pres@gsLineDashPattern  = 1
      ;   mkres@gsMarkerColor = col(yr-years(0)-10)
       ;else if (yr.ge.2000.and.yr.lt.2010) then
        ; pres@gsLineColor  = col(yr-years(0)-20)         ; color of lines"
         ;pres@gsLineDashPattern  = 2
      ;   mkres@gsMarkerColor = col(yr-years(0)-20)
       ;else if (yr.ge.2010) then
      ;   pres@gsLineColor  = col(yr-years(0)-30)         ; color of lines"
      ;   pres@gsLineDashPattern  = 4
      ;   mkres@gsMarkerColor = col(yr-years(0)-30)
      ; else
      ;   res@gsLineColor  = col(yr-years(0))         ; color of lines"
      ;   pres@gsLineDashPattern  = 0
      ;   mkres@gsMarkerColor = col(yr-years(0))
      ; end if
      ; end if
      ; end if
         pres@gsLineColor = col(yr-years(0))         ; color of lines"
         pres@gsLineDashPattern  = 0
         mkres@gsMarkerColor = col(yr-years(0))



       idum=idum+1
       dum(idum) = gsn_add_polyline(wks, plot, raylon, raylat, pres)
       ; Add markers to each of the 'npts' on the great circle path.
       ; mkres@gsMarkerColor  = col(yr-years(0))
       mark(idum) = gsn_add_polymarker(wks,plot,raylon, raylat, mkres)

       lres@gsMarkerColor  = col(yr-years(0))
       gsn_polymarker_ndc(wks,xleg(yr-years(0)),yleg(yr-years(0)),lres)
       gsn_text_ndc      (wks,yr,xtxt(yr-years(0)),ytxt(yr-years(0)),txres)



       ;interpolate rays onto finer resolution for mapping
       do n=0,npts-2
       ;print("n="+n)
        if(.not.ismissing(raylon(n+1)))then
         ngc = ngc0
         gcdist = gc_latlon(raylat(n),raylon(n),raylat(n+1),raylon(n+1),ngc,2)

         if(gcdist/ngc.gt.1.) then
           ;print (raylat(n)+"  "+raylon(n)+"  "+raylat(n+1)+"  "+raylon(n+1))
           ;print ("!!!!ERROR: the distance between 2 points of interpolation is too big ")
           ;print (gcdist)
           ngc = 180
           gcdist = gc_latlon(raylat(n),raylon(n),raylat(n+1),raylon(n+1),ngc,2)
           ;status_exit(1)
         end if
         ;print (gcdist)
         ;print (gcdist@gclat+"  "+gcdist@gclon )  ; print the lats/lons

         do i = 0,ngc-1
          if(gcdist@gclon(i).lt.-0.5) then
           gcdist@gclon(i)=gcdist@gclon(i)+360
          else if(gcdist@gclon(i).gt.359.5)then
           gcdist@gclon(i)=gcdist@gclon(i)-360
          end if
          end if
         end do

         pdf = pdfxy(gcdist@gclon, gcdist@gclat, 360, 181, opt)
         ; set coordinates for pdf2
         ;pdf!0 = "lat"
         ;pdf!1 = "lon"

         ;pdf&lat@units="degrees_north"
         ;pdf&lon@units="degrees_east"
         ;plot = gsn_csm_contour_map (wks,pdf, res)
         ;draw(plot)
         ;frame(wks)



          ;do i = 0,359
          ; do j=0,180
          ;  if (pdf(j,i).gt.0) then
          ;   print (i+"   "+j+"   "+pdf(j,i))
          ;   print (pdf@binx_center(i)+"  "+pdf@biny_center(j))
          ;  end if
          ; end do
          ;end do

          rayfreq = where(pdf.gt.0, 1, rayfreq)
          ;do i = 0,359
          ; do j=0,180
          ;  if (pdf(j,i).gt.0) then
          ;   print (i+"   "+j+"   "+rayfreq(j,i))
          ;   print (pdf@binx_center(i)+"  "+pdf@biny_center(j))
          ;  end if
          ; end do
          ;end do
          ;plot = gsn_csm_contour_map (wks,rayfreq, res)
          ;draw(plot)
          ;frame(wks)



        end if
       end do  ; loop through the ray

       ;plot = gsn_csm_contour_map (wks,rayfreq, res)
       ;draw(plot)
       ;frame(wks)


       ;pres@gsLineColor  = col(k-1)         ; color of lines
       ;pres@xyExplicitLegendLabels = (/"line3","line4"/)
       ; dum1 = gsn_add_polyline(wks,plot, gcdist@gclon ,gcdist@gclat ,pres)
       ;dum1(k-1) = gsn_add_polyline(wks, plot, raylon, raylat, pres)

       ; Add markers to each of the 'npts' on the great circle path.
       ;mkres@gsMarkerColor  = col(k-1)
       ; dum2 = gsn_add_polymarker(wks,plot,gcdist@gclon ,gcdist@gclat ,mkres)
       ;dum2(k-1) = gsn_add_polymarker(wks, plot, raylon, raylat, mkres)

       ;lres@gsMarkerColor  = col(k-1)
       ;gsn_polymarker_ndc(wks,xleg(k-1),yleg(k-1),lres)
       ;gsn_text_ndc      (wks,"k"+k,xtxt(k-1),ytxt(k-1),txres)
       ;draw(plot)
       ;frame(wks)
       delete(ray)
       delete(raylon)
       delete(raylat)

       end if  ; if filein exist
      end do  ; root
     end do  ;year

;---Drawing the plot will draw the attached polylines and markers.
  ;plot = gsn_csm_contour_map (wks,rayfreq, res)

  draw(plot)
  frame(wks)

  delete(dum)

; ***  End plotting  ***************************************************************

    end do  ; k
   end do  ;lon
   end do  ;lat
 end do  ;ip
end do   ; ssn

; end do   ; iy

end
