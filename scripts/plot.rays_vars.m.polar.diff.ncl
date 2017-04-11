; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;var = (/"wnd", "hgt", "RWS")
lev = "300"

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
years = new((/2,3,2/),"integer")

ncase = (/3,3/)

;for 72N 40E
years(0,0,:) = (/1981,1985/)
years(0,1,:) = (/2015,2010/)
years(0,2,:) = (/2000,2003/)

;for 72N 160E
years(1,0,:) = (/1987,2014/)
years(1,1,:) = (/1994,1990/)
years(1,2,:) = (/2016,2012/)

;years = (/1980,1980/)

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

res@tiMainFontHeightF  = 0.018
;res@gsnCenterStringFontHeightF = 0.008

res@gsnRightString   = ""

res@mpGreatCircleLinesOn = False

res@cnFillOn   = True	; color fill
res@cnFillPalette = "blwhre" ; set color map

res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
res@cnMinLevelValF       = -30.			; set the minimum contour level
res@cnMaxLevelValF       = 30.			; set the maximum contour level
res@cnLevelSpacingF      = 5.			; set the interval between contours

resH = True
resH = res
resH@cnLevelSelectionMode = "AutomaticLevels"	; manually set the contour levels with the following 3 resources
resH@cnMinLevelValF       = -8200.			; set the minimum contour level
resH@cnMaxLevelValF       = 9800.			; set the maximum contour level
resH@cnLevelSpacingF      = 100.			; set the interval between contours


resRWS = True
resRWS = res
resRWS@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
resRWS@cnMinLevelValF       = -250.			; set the minimum contour level
resRWS@cnMaxLevelValF       = 250.			; set the maximum contour level
resRWS@cnLevelSpacingF      = 50.			; set the interval between contours


;***********************************

resn = True
;resn = res
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
   do ireg = 0,1
   ;do lat  = 1,1  ;dimsizes(lats)-1
   ;do lon  = 1,1  ;dimsizes(lons)-1
   do ic = 0, ncase(ireg)-1

   ;wind
   fin0 = addfile("../output/matlab/yearly/vars/vars"+lev+"."+years(ireg,ic,0)+".nc","r")
   fin1 = addfile("../output/matlab/yearly/vars/vars"+lev+"."+years(ireg,ic,1)+".nc","r")



   ; read wind
   u0 = fin0->u
   v0 = fin0->v

   u1 = fin1->u
   v1 = fin1->v

   du = u1 - u0
   dv = v1 - v0

   hgt0 = fin0->hgt
   RWS0 = fin0->RWS

   hgt1 = fin1->hgt
   RWS1 = fin1->RWS

   ; smoothing
   nsmth = 5; number of smothing attempts
   do is = 0,nsmth-1
    RWS0 = smth9(RWS0, 0.5,  0.25, True)
    RWS1 = smth9(RWS1, .5, 0.25, True)
   end do



   dhgt = hgt1 - hgt0
   dRWS = RWS1 - RWS0


   copy_VarCoords(u0,du)
   copy_VarCoords(v0,dv)


   copy_VarCoords(hgt0,dhgt)
   copy_VarCoords(RWS0,dRWS)

   ; end read wind





    do k =  1, 6
     ;fout = "../output/matlab/yearly/combplots/raysP_"+season(ssn)+years(ireg,ic,0)+"."+years(ireg,ic,1)+"_"+lats(ireg)+"N_"+lons(ireg)+"E"+"_period"+p+"_k"+k
     rayfreq = 0
     dumu = new(dimsizes(roots)*dimsizes(years(ireg,ic,:)), graphic)
     dumv = new(dimsizes(roots)*dimsizes(years(ireg,ic,:)), graphic)
     dumdu = new(dimsizes(roots)*dimsizes(years(ireg,ic,:)), graphic)
     dumdv = new(dimsizes(roots)*dimsizes(years(ireg,ic,:)), graphic)
     marku = dumu
     markv = dumv
     markdu = dumdu
     markdv = dumdv

     ; ***  Plotting  *******************************************************************

       ;res@gsnCenterString    = "("+lats(ireg)+"N; "+lons(ireg)+"E)      period = "+p+"days     root = "+root

       do iv = 0,1

       if (iv.eq.0)then
        fout = "../output/matlab/yearly/combplots/raysP_wind."+season(ssn)+years(ireg,ic,0)+"."+years(ireg,ic,1)+"_"+lats(ireg)+"N_"+lons(ireg)+"E"+"_period"+p+"_k"+k
        print (fout)
        wks  = gsn_open_wks("png",fout)

        res@tiMainString       = "Ray Path, u " +  years(ireg,ic,0)
         plot(0) = gsn_csm_contour_map_polar(wks,u0,res)

         res@tiMainString       = "Ray Path, v " + years(ireg,ic,0)
         plot(3) = gsn_csm_contour_map_polar(wks,v0,res)

         res@tiMainString       = "Ray Paths, u " +  years(ireg,ic,1)
         plot(1) = gsn_csm_contour_map_polar(wks,u1,res)

         res@tiMainString       = "Ray Paths, v" +  years(ireg,ic,1)
         plot(4) = gsn_csm_contour_map_polar(wks,v1,res)

         res@tiMainString       = "Ray Paths, diff u " + years(ireg,ic,1) + " - "+ years(ireg,ic,0)
         plot(2) = gsn_csm_contour_map_polar(wks,du,res)

         res@tiMainString       = "Ray Paths, diff v " + years(ireg,ic,1) + " - "+ years(ireg,ic,0)
         plot(5) = gsn_csm_contour_map_polar(wks,dv,res)
      end if
      if (iv.eq.1)then
        fout = "../output/matlab/yearly/combplots/raysP_hgt_"+nsmth+"sm9RWS."+season(ssn)+years(ireg,ic,0)+"."+years(ireg,ic,1)+"_"+lats(ireg)+"N_"+lons(ireg)+"E"+"_period"+p+"_k"+k
        print (fout)
        wks  = gsn_open_wks("png",fout)


        res@tiMainString       = "Ray Path, hgt300 " +  years(ireg,ic,0)
        resH@cnMinLevelValF       = 8200			; set the minimum contour level
        resH@cnMaxLevelValF       = 9800			; set the maximum contour level
        resH@cnLevelSpacingF      = 100			; set the interval between contours
        plot(0) = gsn_csm_contour_map_polar(wks,hgt0,resH)

        res@tiMainString       = "Ray Path, RWS " + years(ireg,ic,0)
        plot(3) = gsn_csm_contour_map_polar(wks,RWS0,resRWS)

        res@tiMainString       = "Ray Paths, hgt300 " +  years(ireg,ic,1)
        plot(1) = gsn_csm_contour_map_polar(wks,hgt1,resH)

        res@tiMainString       = "Ray Paths, RWS" +  years(ireg,ic,1)
        plot(4) = gsn_csm_contour_map_polar(wks,RWS1,resRWS)

        res@tiMainString       = "Ray Paths, diff hgt300 " + years(ireg,ic,1) + " - "+ years(ireg,ic,0)
        resH@cnMinLevelValF       = -300.			; set the minimum contour level
        resH@cnMaxLevelValF       = 300.			; set the maximum contour level
        resH@cnLevelSpacingF      = 50.			; set the interval between contours

        plot(2) = gsn_csm_contour_map_polar(wks,dhgt,resH)

        res@tiMainString       = "Ray Paths, diff RWS " + years(ireg,ic,1) + " - "+ years(ireg,ic,0)
        plot(5) = gsn_csm_contour_map_polar(wks,dRWS,resRWS)
     end if

       ;******************************
       ; overlay the plots
       ;******************************
       ;overlay(plot(0),plot1)
       ;******************************

     ; *** Map plotted ******************************************************************

     idum = -1
     do iyr = 0,1
      yr = years(ireg,ic,iyr)
      do root =  roots(0),roots(dimsizes(roots)-1)

      fin = "../output/matlab/yearly/rays/ray_"+season(ssn)+yr+"_"+lats(ireg)+"N_"+lons(ireg)+"E"+"_period"+p+"_k"+k+"_root"+root
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
   ;resP@gsnPanelMainString = "Ray Paths for " + years(ireg,ic,0)+" and "+years(ireg,ic,1)     ; new resource added in NCL V6.4.0
   resP@txString           = "Ray Paths for " + years(ireg,ic,0)+" and "+years(ireg,ic,1) + "  k = "+k +"  period = " + p + "( "+ lons(ireg)+"E, " +lats(ireg)+"N )"
   gsn_panel(wks,plot,(/2,3/),resP)               ; now draw as one plot

   ;frame(wks)

   end do ; iv
   delete(dumu)
   delete(dumv)
   delete(dumdu)
   delete(dumdv)



; ***  End plotting  ***************************************************************

    end do  ; k
   end do  ;case
   ;end do  ;lon
   ;end do  ;lat
   end do ; reg
 end do  ;ip
end do   ; ssn


end
