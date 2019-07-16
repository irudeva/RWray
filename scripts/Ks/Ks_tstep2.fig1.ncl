begin
;************************************************
; read in netCDF file
;************************************************

 dataset = "erain"
 dout = "../../output/Ks/Ks_map/"

 nvar = (/"Ks"/)
 lev = 300

 plot = new(9,graphic)
 fig=(/"a)","b)","c)","d)","e)","f)","g)","h)","i)"/)+" "
 clev_u = (/-80., 80., 10./)
 clev_vu = (/0, 2, .2/)
 clev_Ks= (/0., 20., 2./)



do iv = 0,0
yr=1979

;************************************************
; create default plot
;************************************************
wks_type = "png"
wks_type@wkWidth = 2500
wks_type@wkHeight = 2500

wks = gsn_open_wks(wks_type,dout+nvar(iv)+".v2."+dataset+"."+yr)                  ; send graphics to PNG file
res = True
res@gsnDraw            = False
res@gsnFrame           = False

res@mpMinLatF            = -90      ; range to zoom in on
res@mpMaxLatF            =  90.
res@mpMinLonF            =  0.
res@mpMaxLonF            = 360.
res@mpCenterLonF         = 180.

res@mpProjection = "WinkelTripel"
res@cnFillOn             = True

res@txFontHeightF   = 22

res@lbLabelBarOn               = False

; Set resources for when we turn on labelbar later
res@lbTitleFontHeightF         = .03                 ; make title smaller
;res@lbTitleString              = "mm/day"            ; title string
res@lbTitlePosition            = "Bottom"            ; title position
res@lbTitleDirection           = "Across"            ; title direction
res@lbLabelStride              = 2
res@lbLabelFontHeightF         = 0.02
res@pmLabelBarHeightF          = 0.1
;res@lbBoxEndCapStyle           = "TriangleBothEnds"
res@lbOrientation              = "Horizontal"        ; orientation
;************************************************

;************************************************
; data upload
;************************************************

do n = 0,2
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
 a = addfile("../../output/Ks/Ks_map/Ks."+dataset+"."+chtime+"."+yr+".nc","r")

end if

if(n.eq.1)then
; month - Jan
 chtime = "Jan"
 cht1 = yr+"010100"
 cht2  = yr+"013118"
 ; a = addfile("../../output/Ks/Ks_map/Ks."+chtime+yr+"-"+yr+"."+dataset+".nc","r")
 a = addfile("../../output/Ks/Ks_map/Ks."+dataset+"."+chtime+"."+yr+".nc","r")
end if

if(n.eq.0)then
; ssn - JFM
 chtime = "JFM"
 cht1 = yr+"010100"
 cht2  = yr+"033118"
 a = addfile("../../output/Ks/Ks_map/Ks."+dataset+"."+chtime+"."+yr+".nc","r")
end if

; reading Ks

level = a->level
il = ind(level.eq.lev)

;timeKs = a->time
;utc_timeKs = cd_calendar(timeKs, -3)

;t1Ks = ind(utc_timeKs.eq.cht1)
;t2Ks = ind(utc_timeKs.eq.cht2)

;var = dim_avg_n_Wrap(a->$nvar(iv)$(t1Ks:t2Ks,il,:,:),0)

var = dim_avg_n_Wrap(a->$nvar(iv)$(:,il,:,:),0)
;var = a->$nvar(iv)$

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

; plotting
  if(n.eq.2)then
    res@lbLabelBarOn          = True
  end if
  res@gsnRightString  = chtime+yr+"    "

  cmap = read_colormap_file("BlGrYeOrReVi200")
  res@cnFillPalette = cmap(1:100,:)

  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = clev_Ks(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_Ks(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_Ks(2)			; set the interval between contours

  res@lbLabelStride         = 2

  ;print(nvar(iv)+"")
  ;printVarSummary(var)

  print("plotting")
   res@gsnLeftString   = fig(n*3+2)+"Ks"

  ;plot = gsn_csm_contour_map(wks,var(tstp,il,:,:),res)
  plot(n*3+2) = gsn_csm_contour_map(wks,var,res)
  ;---add zonal average
  ; zres  = True
  ; zres@xyLineThicknesses = 2.
  ; zonal = gsn_csm_attach_zonal_means(wks,plot(2),var,zres)

   delete(cmap)
   cmap = read_colormap_file("BlWhRe")
   res@cnFillPalette = cmap(1:100,:)

   res@gsnLeftString   = fig(n*3)+"U (m/s)"

  ;res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnMinLevelValF  = clev_u(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_u(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_u(2)			; set the interval between contours
  plot(n*3) = gsn_csm_contour_map(wks,u,res)

  res@gsnLeftString   = fig(n*3+1)+"V/U"
  v_u = v ; trick
  v_u = where(u.ge.1.,abs(v/u),v@_FillValue)

  delete(cmap)
  delete(res@cnFillPalette)
  cmap = read_colormap_file("CBR_coldhot")
  res@cnFillPalette = cmap
  res@cnMinLevelValF  = clev_vu(0)			; set the minimum contour level
  res@cnMaxLevelValF  = clev_vu(1)			; set the maximum contour level
  res@cnLevelSpacingF = clev_vu(2)			; set the interval between contours

  plot(n*3+1) = gsn_csm_contour_map(wks,v_u,res)
  print("end plotting")
  delete(cmap)
  delete(res@cnFillPalette)


end do  ; ssn

;************************************************
; create panel
;************************************************
  resP                    = True                 ; modify the panel plot
  resP@gsnPanelMainString = "A common title"     ; set main title

  resP@gsnPanelLeft    = 0.02
  resP@gsnPanelRight    = 0.75                 ; draw from left edge to center

  gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot

end do  ; ivar
end
