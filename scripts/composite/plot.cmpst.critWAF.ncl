begin
chmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
ga=9.80665 ;  Gravitational acceleration


;  Dispersion Statistics
;(10)	 HiQuart
; HighOct="+STATX(12)     ; upper octile
nstat = 12


syear = 1979
lyear = 2017

lat1 = (/60,60,60,60,60/)
lat2 = (/80,80,80,80,80/)

lon1 = (/20,20,60,160,300/)
lon2 = (/80,60,100,200,340/)

do ireg = 1,4

lev = 300

dirPlev = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"

;************************************************
; create points for reg box
;************************************************
 ypts = (/ lat1(ireg), lat1(ireg), lat2(ireg), lat2(ireg), lat1(ireg)/)
 xpts = (/ lon1(ireg), lon2(ireg), lon2(ireg), lon1(ireg), lon1(ireg)/)
 ;************************************************

a = addfile("WAF.mon"+syear+"-"+lyear+".nc","r")

date = a->time  ; time@calendar = "standard"
fdate = cd_calendar(date,-1)
; print(fdate)

lon = a->lon
lat = a->lat
level = a->level

ilev  = ind(level.eq.lev)
; ilon1 = ind(lon.eq.lon1(ireg))
; ilon2 = ind(lon.eq.lon2(ireg))
; ilat1 = ind(lon.eq.lat1(ireg))
; ilat2 = ind(lon.eq.lat2(ireg))

fx = a->Fx(:,ilev,:,:)
fy = a->Fy(:,ilev,:,:)
uclim = a->cluvar(:,ilev,:,:)
vclim = a->clvvar(:,ilev,:,:)

zdev = a->Zdev(:,ilev,:,:)  ; zavar = zvar - czvar
; psidev = a3->psidev(17,3,:,:)

cumag = sqrt(fx^2 + fy^2)
cumag@_FillValue = fx@_FillValue
cumag = where(ismissing(fx), cumag@_FillValue, cumag)
copy_VarMeta(fx,cumag)
; Fx = mask(fx, cumag.lt.10,False)
; Fy = mask(fy, cumag.lt.10,False)
; Fx = mask(fx, lon.ge.20.and.lon.lt.60.and.lat.60.and.lat.le.80,False)
Fx = fx(:,{lat1(ireg):lat2(ireg)},{lon1(ireg):lon2(ireg)})
Fy = fy(:,{lat1(ireg):lat2(ireg)},{lon1(ireg):lon2(ireg)})
WAF = cumag(:,{lat1(ireg):lat2(ireg)},{lon1(ireg):lon2(ireg)})
; WAF = cumag
Lon = lon({lon1(ireg):lon2(ireg)})
Lat = lat({lat1(ireg):lat2(ireg)})

;------------------------------------------------
; PDF
;------------------------------------------------
idate = cd_calendar(date,0)

do imon =0,11

mWAF = WAF(ind(idate(:,1).eq.imon+1),:,:)
mdate  = date (ind(idate(:,1).eq.imon+1))
mfdate = fdate(ind(idate(:,1).eq.imon+1))

sumWAF = dim_avg_n_Wrap(mWAF, (/1,2/))
maxWAF = dim_max_n_Wrap(mWAF, (/1,2/))

WAFpdf  = pdfx(    mWAF, 0, False)    ; default number of bins
sumWAFpdf  = pdfx( sumWAF, 0, False)
maxWAFpdf  = pdfx( maxWAF, 0, False)

;------------------------------------------------
; Percentile
;------------------------------------------------

opt = True
opt@PrintStat = True
stata = stat_dispersion(maxWAF, opt )
print (stata(13))  ; hiDec

statb= stat_dispersion(sumWAF, opt )
print (statb(13))  ; hiDec

printVarSummary(sumWAF)

print(sumWAF(ind(sumWAF.ge.statb(nstat))))
print(mfdate(ind(sumWAF.ge.statb(nstat))))



;------------------------------------------------
; Plot each PDF
;------------------------------------------------
nVar    = 6
plot    = new ( nVar, "graphic")
plota   = new ( nVar, "graphic")
dum     = new ( (/4,3/),graphic)


fout = "cmpst_hioct."+chmon(imon)+"."+lon1(ireg)+"_"+lon2(ireg)+"E."+lat1(ireg)+"_"+lat2(ireg)+"N"
wks  = gsn_open_wks ("png",fout)             ; send graphics to PNG file

res  = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineThicknessF       = 2
res@tiYAxisString          = "PDF (%)"

; res@gsnCenterString        = "WAF"
; plot(0) = gsn_csm_xy (wks, WAFpdf@bin_center, WAFpdf, res)
res@gsnLeftString        = "aver WAF"
res@gsnRightString         = "crit WAF (HighQuartile)= "+statb(nstat) ; stat(10) <-> HighQuartile
plot(0) = gsn_csm_xy (wks, sumWAFpdf@bin_center, sumWAFpdf, res)
; res@gsnCenterString        = "maxWAF"
; plot(2) = gsn_csm_xy (wks, maxWAFpdf@bin_center, maxWAFpdf, res)

; resP    = True
; resP@txString        = "PDFs of Three Variables"
; resP@gsnPanelRowSpec = True                   ; tell panel what order to plt
; gsn_panel(wks,plot,(/1,2/),resP)
; frame(wks)
;------------------------------------------------
; Climatology
;------------------------------------------------
nclyrs =  lyear - syear +1
do iclyr = syear,lyear
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

  zlon = dimsizes(zvar&longitude)
  zlat = dimsizes(zvar&latitude)

  if (iclyr.eq.syear)then
    zarr    = new ((/nclyrs,zlat,zlon/),float)
    t2marr  = new ((/nclyrs,zlat,zlon/),float)
    mslparr = new ((/nclyrs,zlat,zlon/),float)
  end if

 print (iclyr-syear +" "+nclyrs)
  zarr(iclyr-syear,:,:)    = zvar(imon,ilev,:,:)
  t2marr(iclyr-syear,:,:)  = tvar(imon,:,:)
  mslparr(iclyr-syear,:,:) = mslp(imon,:,:)

  delete(zvar)
  delete(tvar)
  delete(mslp)
end do

zclim = dim_avg_n_Wrap(zarr, 0)
t2mclim = dim_avg_n_Wrap(t2marr, 0)
mslpclim = dim_avg_n_Wrap(mslparr, 0)

delete(zarr)
delete(t2marr)
delete(mslparr)

;------------------------------------------------
; Composites
;------------------------------------------------
imdate = cd_calendar(mdate,0)
yrs =  imdate(ind(sumWAF.ge.statb(nstat)),0)
nyrs = num(yrs)

do iyr = 0, nyrs-1
  ;  Surface data
  fname = systemfunc("ls "+dirSrfc+"erain.mslp_air.monmean."+yrs(iyr)+".nc")
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
  mslp = mslp/100

  ;  Plev data
  fname = systemfunc("ls "+dirPlev+"erain.hgt_air_wind.monmean."+yrs(iyr)+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"z").eq."short")then
   zvar = short2flt(clfile->z)
  else
   zvar = clfile->z
  end if
  zvar = zvar/ga

  zlon = dimsizes(zvar&longitude)
  zlat = dimsizes(zvar&latitude)

  if (iyr.eq.0)then
    zarr    = new ((/nyrs,zlat,zlon/),float)
    t2marr  = new ((/nyrs,zlat,zlon/),float)
    mslparr = new ((/nyrs,zlat,zlon/),float)
  end if

  printVarSummary(zvar)
  printVarSummary(zvar(imon,ilev,:,:))
  printVarSummary(zarr(iyr,:,:))
  zarr(iyr,:,:)    = zvar(imon,ilev,:,:)
  t2marr(iyr,:,:)  = tvar(imon,:,:)
  mslparr(iyr,:,:) = mslp(imon,:,:)

end do

zav = dim_avg_n_Wrap(zarr, 0)
t2mav = dim_avg_n_Wrap(t2marr, 0)
mslpav = dim_avg_n_Wrap(mslparr, 0)

;------------------------------------------------
; Plotting Composites
;------------------------------------------------
cmpres                = True
cmpres@gsnDraw        = False
cmpres@gsnFrame       = False
; cmpres@mpOutlineOn    = True	                    ; turn the map outline on
cmpres@gsnMaximize    = True

cmpres@tmXTOn         = False
cmpres@tmYROn         = False

;;set area;;
; cmpres@gsnPolar   = "NH"
cmpres@mpMinLatF                   = 0
cmpres@mpMaxLatF                   = 90
cmpres@mpMinLonF                   = 0
cmpres@mpMaxLonF                   = 360

cmpres@cnFillOn                    = True
cmpres@cnLinesOn                   = False

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

cmap = read_colormap_file("BkBlAqGrYeOrReViWh200")
cmpres@cnFillPalette  = cmap(25:120,:)

cmpres@cnFillOpacityF = 0.6

sres = cmpres
sres@cnFillOn = False                       ; do not color fill
sres@cnLinesOn = True                       ; turn the contour lines on
sres@cnLineLabelsOn = True                  ; turn the line labels on
sres@cnLineLabelInterval = 5
sres@cnLineLabelFontHeightF   = 0.005
sres@tiMainString = ""                      ; do not draw a main title
; delete(sres@cnLevels)                       ; needed as the size of the cnLevels attribute is about to change (next line)
; sres@cnLevels = ispan(-5,35,5)              ; set a different set of contour levels
sres@cnInfoLabelOn = False
; sres@cnInfoLabelOrthogonalPosF = 0.1        ; push the InfoLabel down below the label bar

sres@gsnContourZeroLineThicknessF = 0 	    ; eliminates zero contour
sres@gsnContourNegLineDashPattern = 1       ; sets negative contours to dash pattern 1

; delete(sres@gsnPolar)
delete(sres@cnFillOpacityF)

cmpres@gsnLeftString  = "Hgt";+sprinti("%0.6i", level)
za = zav-zclim
copy_VarMeta(zclim,za)
plot(1) = gsn_csm_contour_map(wks,za,cmpres)
cmpres@gsnLeftString  = "2m temp" ;sprinti("%0.6i", fdate)
t2ma = t2mav-t2mclim
copy_VarMeta(t2mclim,t2ma)
plot(2) = gsn_csm_contour_map(wks,t2ma,cmpres)
cmpres@gsnLeftString  = "Mslp"
mslpa = mslpav-mslpclim
copy_VarMeta(mslpclim,mslpa)
plot(3) = gsn_csm_contour_map(wks,mslpa,cmpres)


plota(1) = gsn_csm_contour(wks,zclim,sres)
plota(2) = gsn_csm_contour(wks,t2mclim,sres)
plota(3) = gsn_csm_contour(wks,mslpclim,sres)

overlay (plot(1), plota(1))
overlay (plot(2), plota(2))
overlay (plot(3), plota(3))


;************************************************
; add the box
;************************************************
resp                  = True                      ; polyline mods desired
resp@gsLineColor      = "red"                     ; color of lines
resp@gsLineThicknessF = 3.0                       ; thickness of lines
; resp@gsLineLabelString= "test"                    ; adds a line label string

; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum(i,0)=gsn_add_polyline(wks,plot(1),xpts(i:i+1),ypts(i:i+1),resp)
    dum(i,1)=gsn_add_polyline(wks,plot(2),xpts(i:i+1),ypts(i:i+1),resp)
    dum(i,2)=gsn_add_polyline(wks,plot(3),xpts(i:i+1),ypts(i:i+1),resp)
  end do

;------------------------------------------------
; Panel Plot
;------------------------------------------------
resP    = True
resP@gsnFrame           = False                ; don't advance panel plot
resP@txString        = chmon(imon)+"    ("+lon1(ireg)+"_"+lon2(ireg)+"E ; "+lat1(ireg)+"_"+lat2(ireg)+"N)"
resP@gsnPanelBottom     = 0.05
; resP@gsnPanelRowSpec = True                   ; tell panel what order to plt
; gsn_panel(wks,plot(1:3),(/2,2/),resP)
gsn_panel(wks,plot(1:3),(/3,1/),resP)


; Draw a text string at the bottom
txres               = True
txres@txFontHeightF = 0.015
gsn_text_ndc(wks,"crit WAF (HighQuartile)= "+statb(nstat),0.5,0.02,txres)
print(sprinti("%6.3i",toint(yrs)))
; gsn_text_ndc(wks,sprinti("%6.3i",toint(yrs)),0.5,0.9,txres)

frame(wks)

;
;
delete(yrs)
delete(mdate)
delete(mfdate)
delete(mWAF)
delete(sumWAF)
delete(maxWAF)
delete(WAFpdf)
delete(sumWAFpdf)
delete(maxWAFpdf)
delete(stata)
delete(statb)
delete(plot)

delete(zarr)
delete(t2marr)
delete(mslparr)

end do  ; mon

end do  ; reg

end
