; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

syr = 1983
eyr = 2017

do jyr = syr,eyr

  yr1 = jyr
  yr2 = yr1
  nyr = yr2-yr1+1

  clmon = (/7/)
  chmon = (/"Jul"/)

  XXwnd = 22  ; critical wind

  level = 300

  fout = "../output/wind/frwnd6h_gt"+XXwnd+"."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"

  i = new((/nyr,3/),"integer")

  do iyr = yr1,yr2

    ffrwnd = addfile("../output/wind/frwnd6h_gt"+XXwnd+"."+chmon(0)+iyr+".erain.nc","r")

    ;---------------------------------------------------------------------------------

    ; read wind fr
    if iyr .eq. yr1 then
      lon = ffrwnd->longitude
      lat = ffrwnd->latitude

      lev = ffrwnd->level
      ilev = ind(lev.eq.level)

      frwnd = new((/nyr,dimsizes(lev),dimsizes(lat),dimsizes(lon)/),integer)
      frwnd_vpos = new((/nyr,dimsizes(lev),dimsizes(lat),dimsizes(lon)/),integer)
      frwnd_vneg = new((/nyr,dimsizes(lev),dimsizes(lat),dimsizes(lon)/),integer)

      frwnd!0 = "time"
      frwnd&time = ispan (yr1,yr2,1)
      frwnd!1 = "level"
      frwnd&level = lev
      frwnd!2 = "latitude"
      frwnd&latitude = lat
      frwnd!3 = "longitude"
      frwnd&longitude = lon

    end if
    tmp = ffrwnd->frwnd

    printVarSummary(tmp)
    printVarSummary(frwnd)
    frwnd(iyr-yr1,:,:,:) = ffrwnd->frwnd
    frwnd_vpos(iyr-yr1,:,:,:) = ffrwnd->frwnd_vpos
    frwnd_vneg(iyr-yr1,:,:,:) = ffrwnd->frwnd_vneg


    ; yrs = ffrwnd->time

  end do

    ; end read wind

  ;---------------------------------------------------------------------------------

  plot = new(2, graphic)


  ; resources

  res                    = True
  res@gsnDraw            = False        ; don't draw plot yet
  res@gsnFrame           = False        ; don't advance frame yet

  ;pres@vpHeightF= 0.4                    ; change aspect ratio of plot
  ;res@vpWidthF = 5.

  ; res@gsnPolar   = "NH"

  res@mpMinLatF          = 20    ; range to zoom in on
  res@mpMaxLatF          = 90

  res@tiMainFontHeightF  = 0.013
  ;res@gsnCenterStringFontHeightF = 0.008

  res@gsnRightString   = ""

  res@mpGreatCircleLinesOn = False

  res@cnFillOn   = True	; color fill
  res@cnFillPalette  = "WhBlGrYeRe"
  ; res@vcLevelPalette = "WhBlGrYeRe"
  ; res@cnFillPalette = "blwhre" ; set color map

  ; res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  ; res@cnMinLevelValF       = -30.			; set the minimum contour level
  ; res@cnMaxLevelValF       = 30.			; set the maximum contour level
  ; res@cnLevelSpacingF      = 5.			; set the interval between contours
  ;

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
  ; pres = True

  ; resources for markers
  ; mkres               = True
  ; mkres@gsMarkerIndex = 17     ; Filled circle
  ; mkres@gsMarkerSizeF = 0.02


  ;legend
  ;
  ; lres               = True
  ; lres@gsMarkerIndex = 17     ; Filled circle
  ; lres@gsMarkerSizeF = 0.03
  ;
  ; txres               = True
  ; txres@txFontHeightF = 0.013
  ;
  ; xleg = (/0.7,0.83 /)   ; Location of
  ; ;;xtxt = (/0.25,0.5/)  ; legend markers
  ; xtxt = xleg + 0.07
  ; yleg = (/.05,0.05/)   ; and text
  ; ytxt = yleg
  ;



  wks  = gsn_open_wks("png",fout)


  res@tiMainString       = "Aver freq wind over " + XXwnd + " m/s, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  resn@gsnRightString   = "level " + lev(ilev)

  ; averaging
  avfrwnd = dim_avg_n_Wrap (frwnd(:,ilev,:,:),0)

  plot(0) = gsn_csm_contour_map(wks,avfrwnd,res)
  draw(plot(0))
  frame(wks)




end do

end
