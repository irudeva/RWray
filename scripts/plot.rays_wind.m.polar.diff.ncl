; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

npts = 21
col = (/"cyan1","green","red","orange","blue","hotpink","olivedrab","cadetblue4","coral4","cyan1"/)
;dum2 = dum1

season=(/"DJF","JJA"/)
;period = (/"-60","60","-30","30","Inf"/)
;period = (/"_inf","50","20"/)
period = (/"Inf"/)

;lats = -ispan(-80, 80, 5)
;lons = ispan(60, 360, 90)
lats = (/72,72/)
lons = (/40,160/)

;years = ispan(iy,iy+9,1)

;for 72N 40E
;years = (/1981,1985/)
;years = (/2015,2010/)
;years = (/2000,2003/)

;for 72N 160E
;years = (/1987,2014/)
;years = (/1990,2015/)
;years = (/1994,1990/)
years = (/2016,2012/)
;years = (/2010,2015/)
;years = (/2010,2015/)
;years = (/2010,2015/)
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

;wind
fu0 = addfile("../output/matlab/yearly1/u."+years(0)+".nc","r")
fv0 = addfile("../output/matlab/yearly1/v."+years(0)+".nc","r")
fu1 = addfile("../output/matlab/yearly1/u."+years(1)+".nc","r")
fv1 = addfile("../output/matlab/yearly1/v."+years(1)+".nc","r")



; read wind
u0 = fu0->u
v0 = fv0->v

u1 = fu1->u
v1 = fv1->v

du = u1 - u0
dv = v1 - v0

copy_VarCoords(u0,du)
copy_VarCoords(v0,dv)

; end read wind

plot = new(6, graphic)


; resources

res                    = True
res@gsnDraw            = False        ; don't draw plot yet
res@gsnFrame           = False        ; don't advance frame yet

;pres@vpHeightF= 0.4                    ; change aspect ratio of plot
;res@vpWidthF = 5.

res@gsnPolar   = "NH"

res@mpMinLatF          = 30    ; range to zoom in on
res@mpMaxLatF          = 90

res@tiMainFontHeightF  = 0.013
;res@gsnCenterStringFontHeightF = 0.008

res@gsnRightString   = ""

res@mpGreatCircleLinesOn = False

res@cnFillOn   = True	; color fill
res@cnFillPalette = "blwhre" ; set color map

res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
res@cnMinLevelValF       = -30.			; set the minimum contour level
res@cnMaxLevelValF       = 30.			; set the maximum contour level
res@cnLevelSpacingF      = 5.			; set the interval between contours

;***********************************

resn = True
resn = res
;resn@gsnTickMarksOn   = False	; no tickmarks
;resn@gsnDraw          = False	; don't draw
;resn@gsnFrame         = False	; don't advance frame
;resn@gsnLeftString    = ""	; no titles
resn@gsnRightString   = ""
;resn@tiXAxisString    = ""
;resn@tiYAxisString    = ""
;resn@cnLineThicknessF = 1.5	; thicker contours
;resn@cnLineLabelsOn   = False	; no line labels

;resn@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;resn@cnMinLevelValF  = -30.			; set the minimum contour level
;resn@cnMaxLevelValF  = 30.			; set the maximum contour level
;resn@cnLevelSpacingF = 5.			; set the interval between contours


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

xleg = (/0.7,0.83 /)   ; Location of
;;xtxt = (/0.25,0.5/)  ; legend markers
xtxt = xleg + 0.07
yleg = (/.05,0.05/)   ; and text
ytxt = yleg






do ssn = 0,0
 do ip    = 0, dimsizes(period)-1
  p = period(ip)
   do lat  = 1,1  ;dimsizes(lats)-1
   do lon  = 1,1  ;dimsizes(lons)-1
    do k =  1, 6
     fout = "../output/matlab/yearly1/raysP_"+season(ssn)+years(0)+"."+years(nyr-1)+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k
     rayfreq = 0
     dumu = new(dimsizes(roots)*dimsizes(years), graphic)
     dumv = new(dimsizes(roots)*dimsizes(years), graphic)
     dumdu = new(dimsizes(roots)*dimsizes(years), graphic)
     dumdv = new(dimsizes(roots)*dimsizes(years), graphic)
     marku = dumu
     markv = dumv
     markdu = dumdu
     markdv = dumdv

     ; ***  Plotting  *******************************************************************
       wks  = gsn_open_wks("png",fout)

       ;res@gsnCenterString    = "("+lats(lat)+"N; "+lons(lon)+"E)      period = "+p+"days     root = "+root

       res@tiMainString       = "Ray Path, u " +  years(0)
       plot(0) = gsn_csm_contour_map_polar(wks,u0,res)

       res@tiMainString       = "Ray Path, v " + years(0)
       plot(3) = gsn_csm_contour_map_polar(wks,v0,res)

       res@tiMainString       = "Ray Paths, u " +  years(1)
       plot(1) = gsn_csm_contour_map_polar(wks,u1,res)

       res@tiMainString       = "Ray Paths, v" +  years(1)
       plot(4) = gsn_csm_contour_map_polar(wks,v1,res)

       res@tiMainString       = "Ray Paths, diff u " + years(1) + " - "+ years(0)
       plot(2) = gsn_csm_contour_map_polar(wks,du,res)

       res@tiMainString       = "Ray Paths, diff v " + years(1) + " - "+ years(0)
       plot(5) = gsn_csm_contour_map_polar(wks,dv,res)
       ;******************************
       ; overlay the plots
       ;******************************
       ;overlay(plot(0),plot1)
       ;******************************

     ; *** Map plotted ******************************************************************

     idum = -1
     do iyr = 0,1
      yr = years(iyr)
      do root =  roots(0),roots(dimsizes(roots)-1)

      fin = "../output/matlab/yearly1/ray_"+season(ssn)+yr+"_"+lats(lat)+"N_"+lons(lon)+"E"+"_period"+p+"_k"+k+"_root"+root
      print(fin)
      if (fileexists(fin)) then
       print("Raypath from "+fin)
       ray = asciiread (fin,(/npts,15/),"float")
       raylon = ray(:,5)
       raylat = ray(:,6)
       ;print(raylon)
       ;print(raylat)


        pres@gsLineColor = col(iyr)         ; color of lines"
        pres@gsLineDashPattern  = 0
        mkres@gsMarkerColor = col(iyr)



       idum=idum+1
       if(iyr.eq.0)then
       print(iyr+" "+idum+" "+root)
        dumu(idum) = gsn_add_polyline(wks, plot(0), raylon, raylat, pres)
        dumv(idum) = gsn_add_polyline(wks, plot(3), raylon, raylat, pres)
        dumdu(idum) = gsn_add_polyline(wks, plot(2), raylon, raylat, pres)
        dumdv(idum) = gsn_add_polyline(wks, plot(5), raylon, raylat, pres)
        ; Add markers to each of the 'npts' on the great circle path.
        ; mkres@gsMarkerColor  = col(yr-years(0))
        marku(idum) = gsn_add_polymarker(wks,plot(0),raylon, raylat, mkres)
        markv(idum) = gsn_add_polymarker(wks,plot(3),raylon, raylat, mkres)
        markdu(idum) = gsn_add_polymarker(wks,plot(2),raylon, raylat, mkres)
        markdv(idum) = gsn_add_polymarker(wks,plot(5),raylon, raylat, mkres)
      end if
      if(iyr.eq.1)then
      print(iyr+" "+idum+" "+root)
       dumu(idum) = gsn_add_polyline(wks, plot(1), raylon, raylat, pres)
       dumv(idum) = gsn_add_polyline(wks, plot(4), raylon, raylat, pres)
       dumdu(idum) = gsn_add_polyline(wks, plot(2), raylon, raylat, pres)
       dumdv(idum) = gsn_add_polyline(wks, plot(5), raylon, raylat, pres)
       ; Add markers to each of the 'npts' on the great circle path.
       ; mkres@gsMarkerColor  = col(yr-years(0))
       marku(idum) = gsn_add_polymarker(wks,plot(1),raylon, raylat, mkres)
       markv(idum) = gsn_add_polymarker(wks,plot(4),raylon, raylat, mkres)
       markdu(idum) = gsn_add_polymarker(wks,plot(2),raylon, raylat, mkres)
       markdv(idum) = gsn_add_polymarker(wks,plot(5),raylon, raylat, mkres)
      end if
       lres@gsMarkerColor  = col(iyr)
       gsn_polymarker_ndc(wks,xleg(iyr),yleg(iyr),lres)
       gsn_text_ndc      (wks,yr,xtxt(iyr),ytxt(iyr),txres)

       delete(ray)
       delete(raylon)
       delete(raylat)

       end if  ; if filein exist
      end do  ; root
     end do  ;year

;---Drawing the plot will draw the attached polylines and markers.
  ;plot = gsn_csm_contour_map (wks,rayfreq, res)

  ;draw(plot(0))


  ;************************************************
  ; create panel
  ;************************************************
   resP                    = True                 ; modify the panel plot
   ;resP@gsnPanelMainString = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
   resP@txString           = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
   gsn_panel(wks,plot,(/2,3/),resP)               ; now draw as one plot

   ;frame(wks)

   delete(dumu)
   delete(dumv)
   delete(dumdu)
   delete(dumdv)



; ***  End plotting  ***************************************************************

    end do  ; k
   end do  ;lon
   end do  ;lat
 end do  ;ip
end do   ; ssn


end
