; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

npts = 20
col = (/"black","red","orange","blue","paleturquoise","green","olivedrab"/)
dum1 = new(6, graphic)
dum2 = dum1

bgs="JJA"

;lon = (/"240","0"/)
;lat = (/"40","30"/)
;print(lon+" "+lat)
;lat = -ispan(-80,80,10)
;lon = lat*0+150
lat = ispan(25,25,10)
lon = ispan(170,170,30)


  do ilat  = 0,dimsizes(lat)-1
  do ilon  = 0,dimsizes(lon)-1

   do root = 3,3
    ;period = (/"-60","60","-30","30","Inf"/)
    ;period = (/"_inf","50","20"/)
    period = (/"Inf"/)
    do ip    = 0, dimsizes(period)-1
     p = period(ip)
     fout = "../output/matlab/poster/wnum_"+bgs+"_"+lat(ilat)+"N_"+lon(ilon)+"E"+"_period"+p+"_root"+root

; ***  Plotting  *******************************************************************
  wks  = gsn_open_wks("png",fout)

  pres = True

   pres@tiMainString       = "k and l wavenumbers"
   pres@tiMainFontHeightF  = 0.013
   ;pres@xyCenterString    = "("+lat(ilat)+"N; "+lon(ilon)+"E)      period = "+p+"days     root = "+root
   ;pres@xyCenterStringFontHeightF = 0.008


  pres@gsnDraw            = False        ; don't draw plot yet
  pres@gsnFrame           = False        ; don't advance frame yet

  pres@vpHeightF= 0.4                    ; change aspect ratio of plot
  pres@vpWidthF = 0.7


  ;pres@trYMinF  = -25.0                   ; min value on y-axis
  ;pres@trYMaxF  =  25.0                   ; max value on y-axis

  pres@tiXAxisFontHeightF = 0.014
  pres@tiYAxisFontHeightF = 0.014
  pres@tmXBLabelFontHeightF = 0.014
  pres@tmYLLabelFontHeightF = 0.014

; legend

   ;lres               = True
   ;lres@xyMarkerIndex = 17     ; Filled circle
   ;lres@xyMarkerSizeF = 0.03

   ;txres               = True
   ;txres@txFontHeightF = 0.013

   xleg = (/0.15,0.15,0.35,0.35,0.56,0.56/)   ; Location of
   xtxt = xleg + 0.07
   yleg = (/0.16,0.10,0.16,0.10,0.16,0.10/)   ; and text
   ytxt = yleg

; *** Map plotted ******************************************************************
  k0 = 0
     do k    = 1,6
     print("k="+k)

      fin = "../output/matlab/ray_"+bgs+"_"+lat(ilat)+"N_"+lon(ilon)+"E"+"_period"+p+"_k"+k+"_root"+root
      ;fin = "../output/matlab/ray_"+bgs+"_lat"+lat(ilat)+"N_lon"+lon(ilon)+"E"+"_period"+p+"_k"+k+"_root"+root
      print(fin)
      if (fileexists(fin)) then
      print("Raypath from "+fin)
       ray = asciiread (fin,(/npts,11/),"float")
       rayk = ray(:,7)
       rayl = ray(:,8)
       rayd = ray(:,2)
       rayd@long_name = "days"
       print(rayd)
       print(rayd)

       if (k0.eq.0) then
         wk = new((/6,dimsizes(rayd)/),float)
         wk@_FillValue = -999
         wk = wk@_FillValue
         wl = wk
         k0 = 1
       end if

       wk(k-1,:)=rayk
       wl(k-1,:)=rayl

       end if  ; if filein exist
      end do  ; k

      pres@xyLineColors  = col
      pres@xyDashPattern = 0
      pres@xyLineThicknesses = 2*(/1.,1,1,1,1,1/)


      pres@tiXAxisString = rayd@long_name
      pres@tiYAxisString = "k"

      pres@pmLegendDisplayMode    = "Always"            ; turn on legend
      pres@pmLegendSide           = "Top"               ; Change location of
      pres@pmLegendParallelPosF   = .2                  ; move units right
      pres@pmLegendOrthogonalPosF = -0.6                 ;
      ;pres@pmLegendOrthogonalPosF = -1.04                 ; move units down
      pres@pmLegendWidthF         = 0.15                ; Change width and
      pres@pmLegendHeightF        = 0.18                ; height of legend.
      pres@lgPerimOn              = False               ; turn off box around
      pres@lgLabelFontHeightF     = .014                 ; label font height
      pres@xyExplicitLegendLabels = "  "+(/"1","2","3","4","5","6"/)         ; create explicit labels


      ;pres@xyExplicitLegendLabels = (/"line3","line4"/)
      ;plot  = gsn_csm_xy (wks,rayd,wk,pres)

      lres = True
       lres = pres
       lres@xyDashPattern = 1
       lres@pmLegendDisplayMode    = "Never"
       lres@tiYAxisString = "l"

       plot  = gsn_csm_xy2 (wks,rayd,wk,wl,pres,lres)



      ;do k    = 1,6
       ;print(k)
       ;lres@gsMarkerColor  = col(k-1)
       ;gsn_polymarker_ndc(wks,xleg(k-1),yleg(k-1),lres)
       ;gsn_text_ndc      (wks,"k"+k,xtxt(k-1),ytxt(k-1),txres)
       ;draw(plot)
      ;end do
      ;frame(wks)
      delete(ray)
      delete(rayk)
      delete(rayd)


;---Drawing the plot will draw the attached polylines and markers.
  draw(plot)
  frame(wks)

; ***  End plotting  ***************************************************************

   end do  ; preiod ip
  end do  ;root
 end do  ;lon
 end do  ;lat

end
