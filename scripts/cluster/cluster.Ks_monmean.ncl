; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

colors = (/"red","green","blue","SkyBlue","yellow"/)
ga=9.80665 ;  Gravitational acceleration

dirPlev = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"
dirKs   = "/Users/irudeva/work/Projects/RWray/output/Ks/Ks_nc/"

yr1 = 1979
yr2 = 2017
nyr = yr2-yr1+1

clrs = new(nyr,string)

clmon = (/8/)
chmon = (/"Aug"/)

XXwnd = 22  ; critical wind

level = 300
nlat = 32

neof =4
ncluster =5
; fout = "../output/wind/frwnd6h_gt"+XXwnd+"."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_clst   = "../output/wind/Ks_p.clst6h_gt"+XXwnd+"."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_eof    = "../output/wind/Ks_p.eof6h_gt"+XXwnd+"."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_z      = "../output/wind/cmpst_p.z_Ksmonmean."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_mslp   = "../output/wind/cmpst_p.mslp_Ksmonmean."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_t2m    = "../output/wind/cmpst_p.t2m_Ksmonmean."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"
fout_frwnd     = "../output/wind/cmpst_p.frwnd_Ksmonmean."+chmon(0)+yr1+"_"+yr2+"."+level+"hPa.erain"

; i = new((/nyr,3/),"integer")

do iyr = yr1,yr2

  ffrwnd = addfile("../output/wind/wind_nc/frwnd6h_gt"+XXwnd+"."+chmon(0)+iyr+".erain.nc","r")

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

  ; printVarSummary(tmp)
  ; printVarSummary(frwnd)
  frwnd(iyr-yr1,:,:,:) = ffrwnd->frwnd
  frwnd_vpos(iyr-yr1,:,:,:) = ffrwnd->frwnd_vpos
  frwnd_vneg(iyr-yr1,:,:,:) = ffrwnd->frwnd_vneg


  ; yrs = ffrwnd->time

end do

  ; end read wind

;------------------------------------------------
; Climatology
;------------------------------------------------
nclyrs =  yr2 - yr1 +1
do iclyr = yr1,yr2
  ;  Surface data
  fname = systemfunc("ls "+dirSrfc+"erain.mslp_air.monmean."+iclyr+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"t2m").eq."short")then
   tvar = short2flt(clfile->t2m)
  else
   tvar = clfile->t2m
  end if
  tvar = tvar + 273.15

  if (getfilevartypes(clfile,"msl").eq."short")then
   mslp = short2flt(clfile->msl)
  else
   mslp = clfile->msl
  end if
  mslp = mslp/100 ; Pa -> hPa

  ;  Plev data
  fname = systemfunc("ls "+dirPlev+"erain.hgt_air_wind.monmean."+iclyr+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"z").eq."short")then
   zvar = short2flt(clfile->z)
  else
   zvar = clfile->z
  end if
  zvar = zvar/ga

  ; Ks data
  ; fname = systemfunc("ls "+dirKs+"erain.Ks_wind.6h."+iclyr+".nc")
  fname = systemfunc("ls "+dirKs+"erain.Ks_wind.monmean."+iclyr+".nc")
  print ("reading Ks file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"Ks").eq."short")then
   ; Ks6h = short2flt(clfile->Ks)
   Ksmon_raw = short2flt(clfile->Ks)
  else
   ; Ks6h = clfile->Ks
   Ksmon_raw = clfile->Ks
  end if
  Ks_lev = clfile->level
  iKslev = ind(doubletofloat(Ks_lev).eq.level)
  ; Ksmon = calculate_monthly_values (Ks6h(:,iKslev,:,:), "avg", 0,False)
  Ksmon = Ksmon_raw(:,iKslev,:,:)

  zlon = dimsizes(zvar&longitude)
  zlat = dimsizes(zvar&latitude)
  zlev = clfile->level
  icllev = ind(lev.eq.level)


  if (iclyr.eq.yr1)then
    zarr    = new ((/nclyrs,zlat,zlon/),float)
    t2marr  = new ((/nclyrs,zlat,zlon/),float)
    mslparr = new ((/nclyrs,zlat,zlon/),float)
    Ksarr   = new ((/nclyrs,zlat,zlon/),float)
    Ksarr1   = new ((/nclyrs,zlat,zlon/),double)
  end if

 print (iclyr-yr1 +" "+nclyrs)
 printVarSummary(zvar)
  imon = clmon(0) -1
  zarr(iclyr-yr1,:,:)    = zvar(imon,icllev,:,:)
  t2marr(iclyr-yr1,:,:)  = tvar(imon,:,:)
  mslparr(iclyr-yr1,:,:) = mslp(imon,:,:)

  Ksarr(iclyr-yr1,:,:)   = doubletofloat(Ksmon(imon,:,:))
  ; to get VarMeta
  if (iclyr.eq.yr1)then
    Ksarr1(iclyr-yr1,:,:) = Ksmon(imon,:,:)
    copy_VarMeta(Ksarr1,Ksarr)
    delete (Ksarr1)
  end if

  delete(zvar)
  delete(tvar)
  delete(mslp)
  ; delete(Ks6h)
  delete(Ksmon)
  delete(Ks_lev)
end do

; copy_VarMeta(Ksmon)
zclim = dim_avg_n_Wrap(zarr, 0)
t2mclim = dim_avg_n_Wrap(t2marr, 0)
mslpclim = dim_avg_n_Wrap(mslparr, 0)
Ksclim = dim_avg_n_Wrap(Ksarr, 0)

; delete(zarr)
; delete(t2marr)
; delete(mslparr)

;------clustering and eof---------------------------------------------------------------
; printVarSummary(frwnd)
Ks1 = Ksarr(lat|:nlat,lon|:,time|:)

clcntr = kmeans_as136(Ks1, ncluster, False)  ; use default options

clcntr!1 = "latitude"
clcntr&latitude = lat(:nlat)
clcntr!2 = "longitude"
clcntr&longitude = lon

printVarSummary (clcntr)

ev     = eofunc(Ks1,neof,False)
ev!1 = "latitude"
ev&latitude = lat(:nlat)
ev!2 = "longitude"
ev&longitude = lon

Ksav      = new ((/ncluster,zlat,zlon/),float)
frwndav      = new ((/ncluster,zlat,zlon/),float)

;---------------------------------------------------------------------------------

plot = new(ncluster+2, graphic)
; resources

res                    = True
res@gsnDraw            = False        ; don't draw plot yet
res@gsnFrame           = False        ; don't advance frame yet

;pres@vpHeightF= 0.4                    ; change aspect ratio of plot
;res@vpWidthF = 5.

res@gsnPolar   = "NH"

res@mpMinLatF          = 20    ; range to zoom in on
res@mpMaxLatF          = 90

res@tiMainFontHeightF  = 0.013
;res@gsnCenterStringFontHeightF = 0.008

res@gsnRightString   = ""

res@mpGreatCircleLinesOn = False

res@cnFillOn   = True	; color fill
res@cnFillPalette  = "testcmp"
; res@cnFillPalette  = "WhBlGrYeRe"
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

;start comment
wks  = gsn_open_wks("PS",fout_clst)


resn@gsnRightString   = "level " + lev(ilev)

;climatology of wind freq
gsn_define_colormap(wks,"testcmap")

do ic = 0,ncluster-1
  icluster = ic +1
  res@tiMainString       = "Cluster "+icluster+", Ks monmean, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"

  ; plot(ic) = gsn_csm_contour_map_polar(wks,clcntr(ic,:,:),res)
  print(clcntr@id)
  ysel = ind(clcntr@id.eq.ic+1)
  print (ysel)
    status_exit(1)
  Ksav(ic,:,:)   = dim_avg_n_Wrap(Ksarr(ysel,:,:), 0)
  Ksa = Ksav(ic,:,:)-Ksclim
  copy_VarMeta(Ksclim,Ksa)
  plot(ic) = gsn_csm_contour_map_polar(wks,Ksa,res)
  delete(ysel)

end do

bres = True

bres@gsnDraw            = False        ; don't draw plot yet
bres@gsnFrame           = False        ; don't advance frame yet


bres@vpHeightF  = 0.18

bres@trYMinF   = 0       ; min value on y-axis
bres@trYMaxF   =  1.5       ; max value on y-axis

bres@tiYAxisString = "cluster" ; y-axis label
bres@tiXAxisString = ""
; bres@tiMainString  = "Darwin Southern Oscillation Index" ; title

bres@gsnYRefLine           = 0.              ; reference line
bres@gsnXYBarChart         = True            ; create bar chart


bars = clcntr@id
print (ispan(yr1,yr2,1))
print (bars)
do j = 0,nyr-1
  clrs(j) = colors(bars(j)-1)
end do
bres@gsnXYBarChartColors = clrs

bres@gsnXYBarChartBarWidth = 0.25

plot(ncluster) = gsn_csm_xy (wks,ispan(yr1,yr2,1),.6 + bars*.2,bres)

lbres = True
lbres@lbPerimOn            = False               ; no label bar box
 ; lbres@lbOrientation        = "Horizontal"        ; orientation
 lbres@vpWidthF             = 0.2                ; size
 lbres@vpHeightF            = 0.1
 lbres@lbLabelFontHeightF   = 0.012               ; label font height
 lbres@lbLabelAlignment     = "InteriorEdges"     ; where to label
 lbres@lbBoxMajorExtentF  = 0.75          ; puts space between color boxes

 lbres@lbMonoFillPattern    = True                ; fill sold
 lbres@lbFillColors         = colors      ; must be RGB triplets
 lbres@lbLabelFontHeightF = 0.008         ; font height. default is small

 gsn_labelbar_ndc (wks,4,(/1,2,3,4/),0.7,0.25,lbres)



;************************************************
; create panel
;************************************************
resP                    = True                 ; modify the panel plot
;resP@gsnPanelMainString = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
; resP@txString           = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
resP@gsnPanelBottom  = 0.05
; gsn_panel(wks,plot,(/ncluster+1,1/),resP)               ; now draw as one plot
gsn_panel(wks,plot(:ncluster+1),(/1+ncluster/2,2/),resP)

yrclust = clcntr@id

print (clcntr@id)
print (clcntr@npts)
print (clcntr@ss2)

status_exit(1)

; frame(wks)
;-----plot eofs------------------------------------------------------------------
print("start eof")
wks1  = gsn_open_wks("PS",fout_eof)


resn@gsnRightString   = "level " + lev(ilev)

do iev = 0,neof-1
  nev = iev+1
  res@tiMainString       = "EOF "+nev+" freq wind over " + XXwnd + " m/s, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  resn@gsnLeftString   =  ev@pcvar(iev)
  plot(iev) = gsn_csm_contour_map_polar(wks1,ev(iev,:,:),res)
end do
;************************************************
; create panel
;************************************************
resP                    = True                 ; modify the panel plot
;resP@gsnPanelMainString = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
; resP@txString           = "Ray Paths for " + years(0)+" and "+years(1)     ; new resource added in NCL V6.4.0
; gsn_panel(wks1,plot(:neof-1),(/neof,1/),resP)               ; now draw as one plot
gsn_panel(wks1,plot(:neof-1),(/neof/2,2/),resP)

print (ev@pcvar)

; draw(plot(0))
; frame(wks)

;------------------------------------------------
; Building Composites
;------------------------------------------------
zav      = new ((/ncluster,zlat,zlon/),float)
t2mav    = new ((/ncluster,zlat,zlon/),float)
mslpav   = new ((/ncluster,zlat,zlon/),float)
Ksav     = new ((/ncluster,zlat,zlon/),float)

zstd     = new ((/ncluster,zlat,zlon/),float)
t2mstd   = new ((/ncluster,zlat,zlon/),float)
mslpstd  = new ((/ncluster,zlat,zlon/),float)
Ksstd    = new ((/ncluster,zlat,zlon/),float)

do nc = 1, ncluster
  ysel = ind(yrclust.eq.nc)
  print(nc)
  print (ysel)
  printVarSummary(dim_avg_n_Wrap(zarr(ysel,:,:), 0))
  printVarSummary(zav(nc-1,:,:) )
  zav(nc-1,:,:)   = dim_avg_n_Wrap(zarr(ysel,:,:), 0)
  t2mav(nc-1,:,:) = dim_avg_n_Wrap(t2marr(ysel,:,:), 0)
  mslpav(nc-1,:,:)= dim_avg_n_Wrap(mslparr(ysel,:,:), 0)
  Ksav(nc-1,:,:)  = dim_avg_n_Wrap(Ksarr(ysel,:,:), 0)

  zstd(nc-1,:,:) = dim_stddev_n_Wrap(zarr(ysel,:,:), 0)
  t2mstd(nc-1,:,:) = dim_stddev_n_Wrap(t2marr(ysel,:,:), 0)
  mslpstd(nc-1,:,:) = dim_stddev_n_Wrap(mslparr(ysel,:,:), 0)
  Ksstd(nc-1,:,:) = dim_stddev_n_Wrap(Ksarr(ysel,:,:), 0)

  delete(ysel)
end do


;------------------------------------------------
; Plotting Composites
;------------------------------------------------
wks2  = gsn_open_wks("PS",fout_z)
wks3  = gsn_open_wks("PS",fout_t2m)
wks4  = gsn_open_wks("PS",fout_mslp)
wks5  = gsn_open_wks("PS",fout_Ks)
gsn_define_colormap(wks5,"testcmap")

cmpres                = True
cmpres                = res

cmpres@pmLabelBarWidthF            = 0.4
cmpres@pmLabelBarHeightF           = 0.05
cmpres@pmLabelBarOrthogonalPosF    = 0.1
cmpres@lbLabelFontHeightF          = 0.006
cmpres@lbLabelAngleF               = 45

; res@cnLevelSelectionMode = "ExplicitLevels" ; use explicit levels
; res@cnLevels = ispan(215,265,5)             ; set the contour levels
; res@cnLineLabelsOn = False                  ; do not use line labels
; res@cnFillOn        = True                  ; color fill
; res@cnLinesOn       = False                 ; do not draw contour lines
; res@cnFillPalette   = "BlueDarkRed18"
; res@tiMainString = "T/U @500hPa"            ; set the main title

cmpres@gsnLeftStringFontHeightF = 0.015
cmpres@gsnLeftString  = "" ;sprinti("%0.6i", fdate)
cmpres@gsnRightString = ""                     ; do not draw the right sub-title (otherwise automatically set to temp@units)

delete(cmpres@cnFillPalette)
; cmap = read_colormap_file("BkBlAqGrYeOrReViWh200")
; cmpres@cnFillPalette  = cmap(25:120,:)
colmap = read_colormap_file("testcmap")
cmpres@cnFillPalette  = colmap(35:180,:)

cmpres@cnFillOpacityF = 0.6

sres = True
sres@cnLineLabelsOn       = True
; sres@cnFillOn = False                       ; do not color fill
; sres@cnLinesOn = True                       ; turn the contour lines on
; sres@cnLineLabelsOn = True                  ; turn the line labels on
; sres@cnLineLabelInterval = 5
; sres@cnLineLabelFontHeightF   = 0.005
; sres@tiMainString = ""                      ; do not draw a main title
; ; delete(sres@cnLevels)                       ; needed as the size of the cnLevels attribute is about to change (next line)
; ; sres@cnLevels = ispan(-5,35,5)              ; set a different set of contour levels
; sres@cnInfoLabelOn = False
; ; sres@cnInfoLabelOrthogonalPosF = 0.1        ; push the InfoLabel down below the label bar
;
; sres@gsnContourZeroLineThicknessF = 0 	    ; eliminates zero contour
; sres@gsnContourNegLineDashPattern = 1       ; sets negative contours to dash pattern 1

; delete(sres@gsnPolar)
; delete(sres@cnFillOpacityF)

cmpres@gsnLeftString  = "Hgt";+sprinti("%0.6i", level)
cmpres@tiMainString   = "average Hgt, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
plot(0) = gsn_csm_contour_map_polar(wks2,zclim,cmpres)
do nc=1,ncluster
  cmpres@tiMainString       = "Hgt composite, wind cluster " + nc + ", "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  za = zav(nc-1,:,:)-zclim
  copy_VarMeta(zclim,za)
  copy_VarMeta(zav,zstd)
;  plot(nc) = gsn_csm_contour_map(wks2,za,cmpres)
  plot(nc+1) = gsn_csm_contour_map_overlay(wks2,za,zstd(nc-1,:,:),cmpres,sres)
end do
resP@txString        = ""
resP@gsnPanelBottom  = 0.05
resP@gsnPanelTop  = 0.95
; gsn_panel(wks2,plot(:ncluster),(/ncluster+1,1/),resP)
gsn_panel(wks2,plot(:ncluster+1),(/1+ncluster/2,2/),resP)
; frame(wks2)

cmpres@gsnLeftString  = "t2m";+sprinti("%0.6i", level)
cmpres@tiMainString   = "average t2m, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
plot(0) = gsn_csm_contour_map_polar(wks3,t2mclim,cmpres)
do nc=1,ncluster
  cmpres@tiMainString       = "t2m composite, wind cluster " + nc + ", "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  t2ma = t2mav(nc-1,:,:)-t2mclim
  copy_VarMeta(t2mclim,t2ma)
  copy_VarMeta(t2mav,t2mstd)
  ; plot(nc) = gsn_csm_contour_map(wks4,t2ma,cmpres)
  plot(nc+1) = gsn_csm_contour_map_overlay(wks3,t2ma,t2mstd(nc-1,:,:),cmpres,sres)
end do
; gsn_panel(wks3,plot(:ncluster),(/ncluster+1,1/),resP)
gsn_panel(wks3,plot(:ncluster+1),(/1+ncluster/2,2/),resP)
; frame(wks3)

cmpres@gsnLeftString  = "Mslp";+sprinti("%0.6i", level)
cmpres@tiMainString   = "average Mslp, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
plot(0) = gsn_csm_contour_map_polar(wks4,mslpclim,cmpres)
do nc=1,ncluster
  cmpres@tiMainString       = "Mslp composite, wind cluster " + nc + ", "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  mslpa = mslpav(nc-1,:,:)-mslpclim
  copy_VarMeta(mslpclim,mslpa)
  copy_VarMeta(mslpav,mslpstd)
  ; plot(nc) = gsn_csm_contour_map(wks2,mslpa,cmpres)
  plot(nc+1) = gsn_csm_contour_map_overlay(wks4,mslpa,mslpstd(nc-1,:,:),cmpres,sres)
end do
; gsn_panel(wks4,plot(:ncluster),(/ncluster+1,1/),resP)
gsn_panel(wks4,plot(:ncluster+1),(/1+ncluster/2,2/),resP)
; frame(wks4)

cmpres@gsnLeftString  = "Ks";+sprinti("%0.6i", level)
cmpres@tiMainString   = "average Ks, "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"

cmpres@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
cmpres@cnMinLevelValF  = 0.			; set the minimum contour level
cmpres@cnMaxLevelValF  = 20.			; set the maximum contour level
cmpres@cnLevelSpacingF = 1.

plot(0) = gsn_csm_contour_map_polar(wks5,Ksclim,cmpres)
do nc=1,ncluster
  cmpres@tiMainString       = "Ks composite, wind cluster " + nc + ", "+chmon(0)+yr1+"_"+yr2+" at  "+lev(ilev)+"hPa"
  Ksa = Ksav(nc-1,:,:)-Ksclim

  cmpres@cnMinLevelValF  = -2.5			; set the minimum contour level
  cmpres@cnMaxLevelValF  = 2.5			; set the maximum contour level
  cmpres@cnLevelSpacingF = .25

  sres@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  sres@cnMinLevelValF  = 0.			; set the minimum contour level
  sres@cnMaxLevelValF  = 20			; set the maximum contour level
  sres@cnLevelSpacingF = 1

  copy_VarMeta(Ksclim,Ksa)
  copy_VarMeta(Ksav,Ksstd)
  ; plot(nc) = gsn_csm_contour_map(wks2,mslpa,cmpres)
  delete(cmpres@cnFillPalette)
  colmap = read_colormap_file("testcmap")
  cmpres@cnFillPalette  = colmap(35:180,:)
  plot(nc+1) = gsn_csm_contour_map_overlay(wks5,Ksa,Ksstd(nc-1,:,:),cmpres,sres)
end do
; gsn_panel(wks5,plot(:ncluster),(/ncluster+1,1/),resP)
gsn_panel(wks5,plot(:ncluster+1),(/1+ncluster/2,2/),resP)
; draw(plot(0))
; frame(wks5)


end
