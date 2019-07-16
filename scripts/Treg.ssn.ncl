
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

diro = "../output/Treg/Treg.ssn/"

; detrending
; dt = 0, chdt = "" <-> raw data
; dt = 1, chdt = "dt" <-> detrended data
; dt = 0
; if (dt.eq.0) then
;   chdt = "raw"
; else if (dt.eq.1) then
;   chdt = "dt"
; end if
; end if

; number of days for running average
; ndav = 0

; window for temperature max
; window = 24*(ndav-1)/2
; window = 60  ; ndays
; rwindow = (window/2)*24   ;

;number of anomalous seasons
nmax = 8

ssnstr = (/"JFM","AMJ","JAS","OND"/)
do issn = 2,3

syr = 1979
eyr = 2017
yrs = ispan(syr,eyr,1)
nyr = eyr-syr+1

ssn=ssnstr(issn)
if (ssn.eq."JFM") then
   syr = 1980
   mStrt = "0101" ; -month start
   mLast = "0331" ; -month last
   maxndays = 90
end if
if  (ssn.eq."JAS") then
   mStrt = "0701" ; -month start
   mLast = "0930" ; -month last
   maxndays = 92
end if
if (ssn.eq."AMJ") then
mStrt = "0401" ; -month start
mLast = "0630" ; -month last
maxndays = 91
end if
if (ssn.eq."OND") then
mStrt = "1001" ; -month start
mLast = "1231" ; -month last
maxndays = 92
end if

; climatology

syr_cl = 1981
eyr_cl = 2010
nyrs_cl = eyr_cl - syr_cl +1

; regions

reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
nreg = dimsizes(reg)
; colors = (/"orange", "red", "blue", "deepskyblue",  \
;              "Khaki", "Khaki4","black", "gray59",\
;                           ; "Khaki", "OliveDrab","black", "gray59",\
;              "chartreuse", "cadetblue4" /)

colors = (/"hotpink", "red", "blue", "deepskyblue",  \
             "Khaki", "orange","black", "gray59",\
             "chartreuse", "cyan" /)



lonreg1 = (/ 5, 30, 70, 100,  80, 100, 235, 260, 30, 160/)
lonreg2 = (/25, 50, 90, 120, 100, 120, 255, 280, 70, 200/)
latreg1 = (/40, 45, 50,  55,  25,  22,  45,  35, 70,  70/)
latreg2 = (/55, 60, 65,  70,  45,  42,  60,  50, 80,  80/)



Tav = new((/nyr,nreg/),float,-999)   ; spacially averaged temperature for each region
; Tav3d = new((/1,nreg,maxndays/),float,-999)   ; spacially averaged temperature for each region
Tav@_FillValue = -999
dtTav = Tav   ; detrended Tav
; time_sel = new((/1,maxndays/),float,-999)
; time_sel@_FillValue = -999

; *********************************************
;  climatology
; *********************************************
print("Calculating climatology")
; do yr = syr_cl,eyr_cl
do yr = syr,eyr
   print(ssn+": reading "+yr)

  ; fin = "/Users/irudeva/work/DATA/ERAint/t2m/eraint.t2m.6h.JAS"+yr+".nc"  ; 6houly
  fin = "../data/t2m/erain.t2m.daily."+ssn+yr+".nc"  ; daily data - after cdo daymean
  a = addfile(fin,"r")

  if (getfilevartypes(a,"t2m").eq."short")then
   tvar = short2flt(a->t2m)
  else
   tvar = a->t2m
  end if
  tvar = tvar - 273.15
  tvar@units = "degC"

  if (getfilevartypes(a,"time").eq."double")then
     time = dble2flt(a->time)
  end if


  ; if (yr.eq.syr)then
  ;  dimt = dimsizes(tvar)
  ;  tvar_all = new((/nyr,dimt(0),dimt(1),dimt(2)/),typeof(tvar))
  ; end if
  ;
  ; tvar_all(yr-syr,:,:,:) = tvar

  ; select season

  ymStrt = yr+mStrt ; year-month start
  ymLast = yr+mLast; year-month last

  YYYYMMDD  = cd_calendar(time, -2) ; convert
  ; YYYYMMDD1 = cd_calendar(time, 0) ; convert
  ; numtimes=dimsizes(time)

  iStrt = ind(ymStrt.eq.YYYYMMDD) ; index of start time
  iLast = ind(ymLast.eq.YYYYMMDD) ; last time


  ; time_sel(0,:dimsizes(time(iStrt:iLast))-1) = time(iStrt:iLast)
  t_tsel = tvar(iStrt:iLast,:,:)

  ; select region
  do ireg = 0, nreg-1
    t_tsel_reg = t_tsel(:,{latreg1(ireg):latreg2(ireg)},{lonreg1(ireg):lonreg2(ireg)})
    ; Tav(yr-syr,ireg,:dimsizes(time(iStrt:iLast))-1) = dim_avg_n(t_rav_rtsel, (/1,2/))

    Tav(yr-syr,ireg) = avg(t_tsel_reg)

    delete(t_tsel_reg)
  end do

delete(tvar)
delete(time)
delete(YYYYMMDD)
delete(t_tsel)
end do

print(Tav(:,0))
;
;
; Tav!0 = "years"
; Tav&years = yrs
;
; tclim = dim_avg_n_Wrap(tvar_all({syr_cl:eyr_cl},:,:,:),0)

; detrending!
; if (dt.eq.1)then
  print("detrending...")
  tavg = dim_avg_n(Tav, 0)
  printVarSummary(tavg)
  dtTav  = dtrend_n(Tav,False,0)
  print(dtTav)
  dtTav  = dtTav + conform(dtTav,tavg,(/1/))
  delete(tavg)
; end if

print("Climatology finished")
; *********************************************
;  read data
; *********************************************

; do yr = syr,eyr
;
;   ; fin = "/Users/irudeva/work/DATA/ERAint/t2m/eraint.t2m.6h.JAS"+yr+".nc"  ; 6houly
;   fin = "../data/t2m/erain.t2m.daily."+ssn+yr+".nc"  ; daily data - after cdo daymean
;   a = addfile(fin,"r")
;   ;
;   ; if (getfilevartypes(a,"t2m").eq."short")then
;   ;  tvar = short2flt(a->t2m)
;   ; else
;   ;  tvar = a->t2m
;   ; end if
;     ; tvar = tvar - 273.15 - tclim ; to degC and minus clim mean
;   tvar = tvar_all(yr-syr,:,:,:) - tclim ; raw/detrended series (if dt==0/1) minus clim mean
;   tvar@units = "degC"          ; update appropriate attribute
;
;   lon  = a->lon
;   lat  = a->lat
;   if (getfilevartypes(a,"time").eq."double")then
;      time = dble2flt(a->time)
;   end if
;
;   ;print(time)
;
;   ; daily values
;   ; opt = True
;   ; opt@nval_crit = 4    ; require at least 4 values per day for a daily mean to be calculated.
;
;   ; tDay = calculate_daily_values (tvar, "avg", 0, opt)
;   ; printVarSummary(tDay)
;
;   ; running average
;   opt = 0 ; set unsmoothed beginning and end pts to x@_FillValue [most common]
;   t_rav = runave_n_Wrap(tvar,ndav,opt,0)  ;
;   ; t_rav = tvar
;
;   ; select season
;
;   ymStrt = yr+mStrt ; year-month start
;   ymLast = yr+mLast; year-month last
;
;   YYYYMMDD  = cd_calendar(time, -2) ; convert
;   ; YYYYMMDD1 = cd_calendar(time, 0) ; convert
;   ; numtimes=dimsizes(time)
;
;   iStrt = ind(ymStrt.eq.YYYYMMDD) ; index of start time
;   iLast = ind(ymLast.eq.YYYYMMDD) ; last time
;
;
;   time_sel(0,:dimsizes(time(iStrt:iLast))-1) = time(iStrt:iLast)
;   t_rav_tsel = t_rav(iStrt:iLast,:,:)
;
;   ; select region
;   do ireg = 0, nreg-1
;     t_rav_rtsel = t_rav_tsel(:,{latreg1(ireg):latreg2(ireg)},{lonreg1(ireg):lonreg2(ireg)})
;     ; Tav(yr-syr,ireg,:dimsizes(time(iStrt:iLast))-1) = dim_avg_n(t_rav_rtsel, (/1,2/))
;     Tav(yr-syr,ireg,:dimsizes(time(iStrt:iLast))-1) = dim_avg_n(t_rav_rtsel, (/0,1,2/))
;
;     delete(t_rav_rtsel)
;   end do
;
;    ; if(yr.eq.syr)then
;    ;
;    ;   Tav3d(0,:,:) =  Tav(yr-syr,:,:)
;    ;   Tall = Tav3d
;    ;   timeall = time_sel
;    ;
;    ; else
;    ;
;    ;   tmp =Tall
;    ;   delete(Tall)
;    ;   Tav3d(0,:,:) =  Tav(yr-syr,:,:)
;    ;   Tall = array_append_record(tmp, Tav3d, 0)
;    ;   delete(tmp)
;    ;
;    ;    tmp = timeall
;    ;    delete(timeall)
;    ;    timeall =  table_attach_rows(tmp, time_sel, 0)
;    ;    delete(tmp)
;    ;
;    ; end if
;
;   delete(tvar)
;   delete(time)
;   delete(t_rav)
;   delete(YYYYMMDD)
;   ; delete(time_sel)
;   delete(t_rav_tsel)
;   Tav3d = Tav@_FillValue
;
;   time_sel = time_sel@_FillValue
;
; end do
;
; status_exit(0)
;
; ; print(cd_calendar(timeall,-2))
;
; ;std
; ; std for a region
; ; Tall_std = new(nreg,float)
; ; do ireg = 0,nreg-1
; ;   Tall_std(ireg) =  dim_stddev_Wrap(ndtooned(Tall(:,ireg,:)))
; ; end do
; ; std for each day and each region
; Tregstd = dim_stddev_n_Wrap(Tav,0)
; do yr = syr,eyr
;    if(yr.eq.syr)then
;
;      Tav3d(0,:,:) =  Tregstd
;      Tstdall = Tav3d
;
;    else
;
;      tmp =Tstdall
;      delete(Tstdall)
;     ;  printVarSummary(Tav3d)
;     ;  printVarSummary(Tregstd)
;      Tav3d(0,:,:) =  Tregstd
;      Tstdall = array_append_record(tmp, Tav3d, 0)
;      delete(tmp)
;
;    end if
; end do
;
; ;*********************************
; ; plot parameters
; ;*********************************
;
;
; tempmax  = new((/nmax,nreg/),float,time_sel@_FillValue)
; timemax  = new((/nmax,nreg/),typeof(timeall),Tall@_FillValue)
; indmax   = new((/nmax,nreg/),integer,-999)
;
; tempmin  = tempmax
; timemin  = timemax
; indmin   = indmax
; ntempmax = tempmax
; ntempmin = tempmax
;
; ; do ireg = nreg-2,nreg-2
; do ireg = 0,nreg-1
;   fout = diro+chdt+"Tregav."+ndav+"dclimanom."+str_squeeze(reg(ireg))+"."+ssn+syr+"_"+eyr
;   wks          = gsn_open_wks ("png",fout)
;
;   res          = True                   ; plot mods desired
;   res@gsnFrame = False                  ; don't advance frame yet
;
;   res@vpHeightF= 0.4                    ; change aspect ratio of plot
;   res@vpWidthF = 0.8
;
;   res@trYMinF  = -15.0                   ; min value on y-axis
;   res@trYMaxF  =  15.0                   ; max value on y-axis
;
;   res@tmXBLabelFontHeightF = 0.012
;   res@tmYLLabelFontHeightF = 0.012
;
;   res@tmXBLabelAngleF      = 90
;
; ; since we stretch the plot, we need to slide it over a bit so it does not
; ; run off the page. we do this by:
;
;   res@vpXF     = 0.1                   ; start plot at x ndc coord
;
;   res@tiMainString     = str_squeeze(reg(ireg))+":  Climaotological anomalies of T, " + ndav + "-day average, " + ssn
;   res@tiMainFontHeightF   = 0.015
;   res@gsnLeftString     = "window = "+window+" days"
;   res@gsnStringFontHeightF   = 0.01
;   res@tiYAxisString    = "T anomalies, ~S~o~N~C"    ; y-axis label
;   res@tiYAxisFontHeightF    = 0.012
;   res@tiXAxisString    = "time"    ; y-axis label
;   res@tiXAxisFontHeightF    = 0.012
;
; ; create a reference line and shade values above and below with
; ; selected colors. This is shading array dsoid.
; ;   res@gsnYRefLine           = 0.0             ; create a reference line
; ;*********************************
; ; polyline parameters used on both plots
; ;*********************************
; ;   polyres                  = True
; ;   polyres@gsLineThicknessF = 3.0
; ;*********************************
; ; first plot in color
; ;*********************************
; ;   res@gsnAboveYRefLineColor = "red"              ; above ref line fill red
; ;   res@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue
;
;   y = ndtooned(Tall(:,ireg,:))
;   ystd = ndtooned(Tstdall(:,ireg,:))
;   ; print(y1(:30))
;   ; print(y(:30))
;   ; status_exit(1)
;   t = ndtooned(timeall)
;   ; printVarSummary(y)
;   x = ispan(1,dimsizes(y),1)
;
;   res@tmXBMode        = "Explicit"
;   res@tmXBValues      = fspan(0,dimsizes(y),nyr+1)
;   ; res@tmXBLabels      = "" + res@tmXBValues
;   res@tmXBLabels      = "" + ispan(syr,eyr+1,1)
;   ; res@tmXBMinorValues = ispan(1949,2006,1)
;   ; res@gsnLeftString   = "tmXBMode = Explicit"
;
;   res@trXMinF  = min(x)-2
;   res@trXMaxF  =  max(x)+2
;
;   plot = gsn_csm_xy (wks,x,y,res)        ; create plot
; ;   gsn_polyline(wks,plot,dateF,(/dsoid/),polyres) ; add polyline
;   frame(wks)                                     ; now advance frame
;
;   ; status_exit(1)
;
; ;*********************************
; ; cut out T max and min
; ;*********************************
;
; ynewmax    = y
; ycutmax = new(dimsizes(y),typeof(y),-999)
; ycutmax = ycutmax@_FillValue
;
; ynewmin = y
; ycutmin = ycutmax
;
; t@calendar =	"standard"
; t@units = "hours since 1900-1-1 00:00:00"
;
; do np = 0,nmax-1
;   ;  printVarSummary(ynew)
;    i = maxind(ynewmax)
;    print(i+"   "+np)
;    print("time of anom: "+t(i) + "  "+cd_calendar(t(i), -2))
;    print("temp anom   = "+ ynewmax(i))
;
;    ynewmax = where(t.gt.t(i)-rwindow .and. t.lt.t(i)+rwindow,Tall@_FillValue,ynewmax)
;    ycutmax = where(t.ge.t(i)-rwindow .and. t.le.t(i)+rwindow,y,ycutmax)
;
;    timemax(np,ireg) = t(i)
;    tempmax(np,ireg) = y(i)
;    ntempmax(np,ireg) = y(i)/ystd(i)
;    indmax(np,ireg)  = i
;
;    j = minind(ynewmin)
;   ;  print(i+"   "+np)
;
;    ynewmin = where(t.gt.t(j)-rwindow .and. t.lt.t(j)+rwindow,Tall@_FillValue,ynewmin)
;    ycutmin = where(t.ge.t(j)-rwindow .and. t.le.t(j)+rwindow,y,ycutmin)
;
;    timemin(np,ireg) = t(j)
;    tempmin(np,ireg) = y(j)
;    ntempmin(np,ireg) = y(j)/ystd(j)
;    indmin(np,ireg)  = j
;
; end do
;
; ; status_exit(1)
;
; pmres = True
; pmres@gsMarkerIndex  = 16
;
; ; print(ycut)
; fout_cut = fout+".cut"
; print(fout_cut)
; wks1   = gsn_open_wks ("png",fout_cut)
; res@xyLineColor  = "blue"
; plot  = gsn_csm_xy (wks1,x,ynewmax ,res)
; ; gsn_polyline(wks1,plot,x,ynewmin,res)
; res1=True
; res1 = res
; res1@gsLineColor   = "red"
; gsn_polyline(wks1,plot,x,ycutmax,res1)
; pmres@gsMarkerColor  = res1@gsLineColor
; gsn_polymarker(wks1,plot,x(indmax(:,ireg)),tempmax(:,ireg),pmres)
;
;
; res1@gsLineColor   = "green"
; gsn_polyline(wks1,plot,x,ycutmin,res1)
; pmres@gsMarkerColor  = res1@gsLineColor
; gsn_polymarker(wks1,plot,x(indmin(:,ireg)),tempmin(:,ireg),pmres)
; frame(wks1)
;
; end do  ; end plotting regional T anomalies
;*********************************
; Write netcdf output
;*********************************
do ireg = 0,nreg-1

  ;filo = "dates.dtTmax."+str_squeeze(reg(ireg))+"."+ssn+syr+"_"+eyr+".nc"             ; Output file - problem with region name
  print(diro)
  system("/bin/rm -f "  + diro + "t2m."+"reg"+ireg+"."+ssn+syr+"_"+eyr+".nc"  )    ; remove if exists

	fout  = addfile(diro+ "t2m."+"reg"+ireg+"."+ssn+syr+"_"+eyr+".nc"   , "c")  ; open output file

  setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
	; fAtt@title         = "NCL Efficient Approach to netCDF Creation"
; 	fAtt@source_file   =  "original-file.nc"
fAtt@creation_date = systemfunc ("date")
; if (dt.eq.0) then
fAtt@Content   = "Seasonal raw and detrended t2m values in the selected region and "+nmax+" years of max/min t2m."
; else if (dt.eq.1) then
  ; fAtt@Content   = "Detrended seasonal t2m values in the selected regionand "+nmax+" years of max/min t2m."
; end if
; end if
fAtt@region = reg(ireg)

print(fAtt)

fileattdef( fout, fAtt )            ; copy file attributes

;===================================================================

dimNames = (/ "year",  "n" /)
dimSizes = (/    nyr, nmax /)
dimUnlim = (/  False, False/)
filedimdef(fout,dimNames,dimSizes,dimUnlim)

;===================================================================
; create variables
;years@region = reg(ireg)
filevardef(fout, "year" ,"short","year" )
filevardef(fout,   "t2m" , "float","year" )
filevardef(fout,   "t2m_dt" , "float","year" )
filevardef(fout,   "yr_t2m_max" , "short","n" )
filevardef(fout,   "yr_t2m_min" , "short","n" )
filevardef(fout,   "yr_t2mdt_max" , "short","n" )
filevardef(fout,   "yr_t2mdt_min" , "short","n" )


;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================


temp = 0.
temp@region = reg(ireg)
temp@units = "deg C"
temp@long_name = ssn+" t2m"
filevarattdef(fout,"t2m",temp)
temp@long_name = "detrended "+ssn+" t2m"
filevarattdef(fout,"t2m_dt",temp)

year = 0
year@region = reg(ireg)
year@long_name = "year"
filevarattdef(fout,"year",year)
year@long_name = "years of max "+ssn+" tempertaures"
filevarattdef(fout,"yr_t2m_max",year)
year@long_name = "years of min "+ssn+" tempertaures"
filevarattdef(fout,"yr_t2m_min",year)
filevarattdef(fout,"year",year)
year@long_name = "years of max "+ssn+" detrended tempertaures"
filevarattdef(fout,"yr_t2mdt_max",year)
year@long_name = "years of min "+ssn+" detrended tempertaures"
filevarattdef(fout,"yr_t2mdt_min",year)


;===================================================================

setfileoption(fout,"DefineMode",False)

;===================================================================

; fout->region   = stringtochar(str_squeeze(reg(ireg))  )
fout->t2m     = (/Tav(:,ireg)/)
fout->t2m_dt     = (/dtTav(:,ireg)/)
fout->year     = (/integertoshort(yrs)/)

indexes = dim_pqsort(Tav(:,ireg), 1)
fout->yr_t2m_max   = (/integertoshort(yrs(indexes(nyr-nmax:)))/)
fout->yr_t2m_min   = (/integertoshort(yrs(indexes(:nmax-1)))/)
indexes = dim_pqsort(dtTav(:,ireg), 1)
fout->yr_t2mdt_max   = (/integertoshort(yrs(indexes(nyr-nmax:)))/)
fout->yr_t2mdt_min   = (/integertoshort(yrs(indexes(:nmax-1)))/)
end do  ; nregs

; ;===================================================================
; ;plot anomalies
; ;===================================================================
;
; dummax=new(nreg-1,graphic)
; dummin=new(nreg-1,graphic)
;
; dummax1=new(nreg-1,graphic)
; dummin1=new(nreg-1,graphic)
;
;
; wks_time   = gsn_open_wks ("png",diro+chdt+"Tregav."+ndav+"dclimanom."+ssn+syr+"_"+eyr)             ; send graphics to PNG file
; print(diro+chdt+"Tregav."+ndav+"dclimanom."+ssn+syr+"_"+eyr)
;
; res@gsnDraw  = False
; res@gsnFrame = False
;
; res@trYMinF  = -20.0                   ; min value on y-axis
; res@trYMaxF  = 20.0                   ; max value on y-axis
;
;
; res1                   = True
;
; if (dt.eq.0)then
;   res@tiMainString      = "Climaotological anomalies of T, " + ndav + "-day average, " + ssn            ; add title
; else if (dt.eq.1)then
;   res@tiMainString      = "Climaotological anomalies of detrended T, " + ndav + "-day average, " + ssn            ; add title
; end if
; end if
; res@xyMarkLineModes   = "Markers"                ; choose which have markers
; res@tmLabelAutoStride = True                     ; nice tick mark labels
;
;
; zmax = ycutmax
; zmax = zmax@_FillValue
; zmax(indmax(:,0)) = tempmax(:,0)
;
; zmin = ycutmin
; zmin = zmin@_FillValue
; zmin(indmin(:,0)) = tempmin(:,0)
;
; res@xyMarker  = 16
; res@xyMarkerColor = colors(0)
; pmres1 = True
; pmres1@gsMarkerIndex  = 16
; pmres1@gsMarkerColor = colors(0)
;
; plot0  = gsn_csm_xy (wks_time,x,zmax,res)
; dum0   = gsn_add_polymarker(wks_time,plot0,x,zmin,pmres1)
;
;
; do ireg = 1,nreg-1
;     zmax = zmax@_FillValue
;     zmax(indmax(:,ireg)) = tempmax(:,ireg)
;     zmin = zmax@_FillValue
;     zmin(indmin(:,ireg)) = tempmin(:,ireg)
;     pmres1@gsMarkerColor  = colors(ireg)
;     pmres1@gsMarkerIndex  = 16
;     dummax(ireg-1)=gsn_add_polymarker(wks_time,plot0,x,zmax,pmres1)
;     pmres1@gsMarkerIndex  = 16
;     dummin(ireg-1)=gsn_add_polymarker(wks_time,plot0,x,zmin,pmres1)
; end do
;
; ; draw(plot0)
;
; ; frame(wks_time)
; ; status_exit(1)
; ;***********************************************
; ; legend resources
; ;***********************************************
;   lgres                    = True
;   lgres@lgMonoItemType     = True
;   lgres@lgMarkerColors     = colors(::-1)
;   lgres@lgItemType         = "Markers"        ; show lines only (default)
;   lgres@lgMonoMarkerIndex  = True
;   lgres@lgMarkerIndex      = 16
;   lgres@lgMonoMarkerSize   = True
;   lgres@lgMarkerSizeF      = .005
;   lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
;   lgres@vpWidthF           = 0.15           ; width of legend (NDC)
;   lgres@vpHeightF          = 0.2            ; height of legend (NDC)
;
;   ; lgres@lgPerimColor       = "orange"       ; draw the box perimeter in orange
;   ; lgres@lgPerimThicknessF  = 5.0            ; thicken the box perimeter
;   lbid = gsn_create_legend(wks_time,nreg,reg(::-1),lgres)         ; create legend
;
;   amres = True
;   amres@amParallelPosF   = 0.66	      	    ; move legend to the right
;   amres@amOrthogonalPosF = -.25           ; move the legend down
;   annoid1 = gsn_add_annotation(plot0,lbid,amres)   ; add legend to plot
;
;
; res1 = res
; res1@trYMinF  = -5.0                   ; min value on y-axis
; res1@trYMaxF  = 5.                   ; max value on y-axis
;
; zmax = zmax@_FillValue
; ; zmax(indmax(:,0)) = tempmax(:,0)/Tall_std(0)
; zmax(indmax(:,0)) = ntempmax(:,0)
; zmin = zmin@_FillValue
; ; zmin(indmin(:,0)) = tempmin(:,0)/Tall_std(0)
; zmin(indmin(:,0)) = ntempmin(:,0)
;
; if (dt.eq.0)then
;   res1@tiMainString     = "Standardized climaotological anomalies of T, " + ndav + "-day average, " + ssn
; else if (dt.eq.1)then
;  res1@tiMainString     = "Standardized climaotological anomalies of detrended T, " + ndav + "-day average, " + ssn
; end if
; end if
; res1@tiYAxisString    = "normalized T anomalies, ~S~o~N~C"
; res1@xyMarkerColor = colors(0)
; plot1 =  gsn_csm_xy (wks_time,x,zmax,res1)
; pmres1@gsMarkerColor = colors(0)
; pmres1@gsMarkerIndex  = 16
; dum01 = gsn_add_polymarker(wks_time,plot1,x,zmin,pmres1)
;
; ; draw(plot1)
; ;  frame(wks_time)
; ; status_exit(1)
;
; do ireg = 1,nreg-1
;     zmax = zmax@_FillValue
;     ; zmax(indmax(:,ireg)) = tempmax(:,ireg)/Tall_std(ireg)
;     zmax(indmax(:,ireg)) = ntempmax(:,ireg)
;     zmin = zmax@_FillValue
;     ; zmin(indmin(:,ireg)) = tempmin(:,ireg)/Tall_std(ireg)
;     zmin(indmin(:,ireg)) = ntempmin(:,ireg)
;
;     pmres1@gsMarkerColor  = colors(ireg)
;     pmres1@gsMarkerIndex  = 16
;     dummax1(ireg-1)=gsn_add_polymarker(wks_time,plot1,x,zmax,pmres1)
;     pmres1@gsMarkerIndex  = 16
;     dummin1(ireg-1)=gsn_add_polymarker(wks_time,plot1,x,zmin,pmres1)
; end do
;
; ;************************************************
; ; create panel
; ;************************************************
; gsn_panel(wks,(/plot0,plot1/),(/2,1/),False)             ; now draw as one plot
;
; frame(wks_time)

end do ; ssn
end
