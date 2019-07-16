begin
;************************************************
; read in netCDF file
;************************************************

 dataset = "erain"
 dout = "../../output/Ks/Ks_map/"

 nvar = (/"Ks"/)
 lev = 300

 plot = new(6,graphic)

do iv = 0,0
yr=1979

do n = 0,3
if(n.eq.3)then
; week - 12-18Jan
 chtime = "12-18Jan"
 cht1 = yr+"011200"
 cht2  = yr+"011818"
end if


if(n.eq.2)then
; day - 15Jan
 chtime = "15Jan"
 cht1 = yr+"011500"
 cht2  = yr+"011518"
end if

if(n.eq.1)then
; month - Jan
 chtime = "Jan"
 cht1 = yr+"010100"
 cht2  = yr+"013118"
end if

if(n.eq.0)then
; ssn - JFM
 chtime = "JFM"
 cht1 = yr+"010100"
 cht2  = yr+"033118"
end if

; ssn = (/"JJA", "DJF"/)

  a = addfile("../../output/Ks/Ks_nc/"+dataset+".Ks_wind.6h."+yr+".nc","r")
  ;************************************************
  if(any((/"Ks"/).eq.nvar(iv))) then
    clev_Ks= (/0., 20., 2./)
  end if
  ;************************************************
  level = a->level
  il = ind(level.eq.lev)

  timeKs = a->time
  utc_timeKs = cd_calendar(timeKs, -3)

  t1Ks = ind(utc_timeKs.eq.cht1)
  t2Ks = ind(utc_timeKs.eq.cht2)

  var = dim_avg_n_Wrap(a->$nvar(iv)$(t1Ks:t2Ks,il,:,:),0)

  if (any(isnan_ieee(var))) then
      value = 1.e20
      replace_ieeenan (var, value, 0)
      var@_FillValue = value
  end if

  ; reading wind

  b=addfile("~/work/DATA/ERAint/Plev/"+dataset+".hgt_air_wind.6h."+yr+".nc","r")
  ulevel = b->level
  ilu = ind(ulevel.eq.lev)

  timeu =b->time
  utc_timeu = cd_calendar(timeu, -3)
  t1u = ind(utc_timeu.eq.cht1)
  t2u = ind(utc_timeu.eq.cht2)

  u = dim_avg_n_Wrap(short2flt(b->u(t1u:t2u,ilu,:,:)),0)
  v = dim_avg_n_Wrap(short2flt(b->v(t1u:t2u,ilu,:,:)),0)
  clev_u = (/-80., 80., 10./)
  clev_vu = (/0, 2, .2/)

  ;************************************************
  ; create default plot
  ;************************************************
  wks_type = "png"
  wks_type@wkWidth = 2500
  wks_type@wkHeight = 2500

  wks = gsn_open_wks(wks_type,dout+nvar(iv)+"."+dataset+"."+chtime+yr)                  ; send graphics to PNG file
  res = True
  res@gsnDraw            = False
  res@gsnFrame           = False

  res@mpMinLatF            = -90      ; range to zoom in on
  res@mpMaxLatF            =  90.
  res@mpMinLonF            =  0.
  res@mpMaxLonF            = 360.
  res@mpCenterLonF         = 180.

  res@mpProjection = "WinkelTripel"

  res@cnFillOn             = True                  ; turn on color fill
  cmap = read_colormap_file("BlGrYeOrReVi200")
  res@cnFillPalette = cmap(1:100,:)

  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = clev_Ks(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_Ks(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_Ks(2)			; set the interval between contours

  ;res@cnSmoothingOn = True

  res@lbLabelStride         = 2

  print(nvar(iv)+"")
  printVarSummary(var)

  print("plotting")
  ; res@gsnZonalMean    = True

  ;plot = gsn_csm_contour_map(wks,var(tstp,il,:,:),res)
  plot(2) = gsn_csm_contour_map(wks,var,res)
  ;---add zonal average
   zres  = True
   zres@xyLineThicknesses = 2.
   zonal = gsn_csm_attach_zonal_means(wks,plot(2),var,zres)

   delete(cmap)
   cmap = read_colormap_file("BlWhRe")
   res@cnFillPalette = cmap(1:100,:)

  ;res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnMinLevelValF  = clev_u(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_u(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_u(2)			; set the interval between contours
  plot(0) = gsn_csm_contour_map(wks,u,res)
  v_u = v ; trick
  v_u = where(u.ge.5.,abs(v/u),v@_FillValue)

  delete(cmap)
  delete(res@cnFillPalette)
  cmap = read_colormap_file("CBR_coldhot")
  res@cnFillPalette = cmap
  res@cnMinLevelValF  = clev_vu(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_vu(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_vu(2)			; set the interval between contours

  plot(1) = gsn_csm_contour_map(wks,v_u,res)
  print("end plotting")
  delete(cmap)
  delete(res@cnFillPalette)



  ;************************************************
  ; create panel
  ;************************************************
    resP                    = True                 ; modify the panel plot
    resP@gsnPanelMainString = "A common title"     ; set main title

    resP@gsnPanelLeft    = 0.02
    resP@gsnPanelRight    = 0.75                 ; draw from left edge to center

    gsn_panel(wks,plot(:2),(/1,3/),resP)               ; now draw as one plot

end do  ; ssn
end do  ; ivar
end
