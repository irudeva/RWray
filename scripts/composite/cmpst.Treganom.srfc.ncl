; function remove29Feb(var,time,yr)
; begin 
;   utc_time = cd_calendar(time, -2)

;   ch28Feb = yr+"0228" 
;   ch1Mar  = yr+"0301"

;   t1 = ind(utc_srfctime.eq.ch28Feb)
;   t2 = ind(utc_srfctime.eq.ch1Mar)

;   timenew = array_append_record(time(:t1),time(t2:),0)
;   varnew = array_append_record(var(:t1,:,:),var(t2:,:,:),0)
;   delete(time)
;   delete(var)

;   time = time
;   var = varnew

;   delete(timenew)
;   delete(varnew)
;   delete(utc_time)

;   return(var, time)
; end


begin

reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
nreg = dimsizes(reg)

ssn = "JAS"
syr0 = 1980
eyr0 = 2017

syr = 1979
eyr = 2017
yrs = ispan(syr,eyr,1)
nyr = eyr-syr+1

syr_cl = 1980 ;1981
eyr_cl = 2001
nyrs_cl = eyr_cl - syr_cl +1

; number of days for running average
ndav = 15

; number of extreme events
nmax = 15

ga=9.80665 ;  Gravitational acceleration

level = 300
;number of days for composits
ndcmpst = 10
tadj1 = 0
tadj2 = ndcmpst+tadj1-1

dirin   = "../output/Treg/" 
diro    = "../output/Treg/composites/"
dirPlev = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"
dirKs   = "/Users/irudeva/work/Projects/RWray/output/Ks/Ks_nc/"

sign  = (/"tmax","tmin"/)
trend = (/"raw","detrended"/)

;------------------------------------------------
; read data
;------------------------------------------------

do iyr = syr,eyr
 
    ch28Feb = iyr+"0228"
    ch1Mar  = iyr+"0301"

  ;  Surface data
  fname = systemfunc("ls "+dirSrfc+"erain.mslp_air.daily."+iyr+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"t2m").eq."short")then
   tvar = short2flt(clfile->t2m)
  else
   tvar = clfile->t2m
  end if
  tvar = tvar - 273.15

  if (getfilevartypes(clfile,"msl").eq."short")then
   mslpvar = short2flt(clfile->msl)
  else
   mslpvar = clfile->msl
  end if
  mslpvar = mslpvar/100. ; Pa -> hPa

  lonsrfc = mslpvar&lon
  latsrfc = mslpvar&lat

  nlonsrfc = dimsizes(mslpvar&lon)
  nlatsrfc = dimsizes(mslpvar&lat)

  ;removing 29Feb
  srfctime = clfile->time
  utc_srfctime = cd_calendar(srfctime, -2)

;   print(cd_calendar(srfctime, -3))

  
  t1 = ind(utc_srfctime.eq.ch28Feb)
  t2 = ind(utc_srfctime.eq.ch1Mar)

  srfctimenew = array_append_record(srfctime(:t1),srfctime(t2:),0)
  mslpvarnew = array_append_record(mslpvar(:t1,:,:),mslpvar(t2:,:,:),0)
  tvarnew = array_append_record(tvar(:t1,:,:),tvar(t2:,:,:),0)
  delete(srfctime)
  delete(mslpvar)
  delete(tvar)

  srfctime = srfctimenew
  mslpvar = mslpvarnew
  tvar = tvarnew

  delete(srfctimenew)
  delete(mslpvarnew)
  delete(tvarnew)
  delete(utc_srfctime)

;*************************************************************************
;  Plev data

;   fname = systemfunc("ls "+dirPlev+"erain.hgt_air_wind.daily."+iyr+".nc")
;   print ("reading file "+ fname)
;   clfile = addfile(fname,"r")

;   zlev = clfile->level
;   icllev = ind(zlev.eq.level)


;   if (getfilevartypes(clfile,"z").eq."short")then
;    zvar = short2flt(clfile->z(:,icllev,:,:))
;   else
;    zvar = clfile->z(:,icllev,:,:)
;   end if
;   zvar = zvar/ga
;   zlon = dimsizes(zvar&lon)
;   zlat = dimsizes(zvar&lat)


;   print("getting wind components")
;   if (getfilevartypes(clfile,"u").eq."short")then
;    uvar = short2flt(clfile->u(:,icllev,:,:))
;   else
;    uvar = clfile->u(:,icllev,:,:)
;   end if
;   u2 = uvar*uvar

;   if (getfilevartypes(clfile,"v").eq."short")then
;    vvar = short2flt(clfile->v(:,icllev,:,:))
;   else
;    vvar = clfile->v(:,icllev,:,:)
;   end if
;   v2 = vvar*vvar
;   print("calculating the wind speed...")
;   wvar = sqrt(u2+v2)
;   delete(uvar)
;   delete(vvar)
;   delete(u2)
;   delete(v2)


;   wvar!0 = "time"
;   wvar&time = clfile->time
;   wvar!1 = "latitude"
;   wvar&latitude = clfile->lat
;   wvar!2 = "longitude"
;   wvar&longitude = clfile->lon
;   print("end calculating the wind speed")

; ;removing 29Feb in Plev
;   ztime = clfile->time
;   utc_ztime = cd_calendar(ztime, -2)
  
;   t1 = ind(utc_ztime.eq.ch28Feb)
;   t2 = ind(utc_ztime.eq.ch1Mar)

;   ztimenew = array_append_record(ztime(:t1),ztime(t2:),0)
;   zvarnew = array_append_record(zvar(:t1,:,:),zvar(t2:,:,:),0)
;   wvarnew = array_append_record(wvar(:t1,:,:),wvar(t2:,:,:),0)
;   delete(ztime)
;   delete(zvar)
;   delete(wvar)

;   ztime = ztimenew
;   zvar = zvarnew
;   wvar = wvarnew

;   delete(ztimenew)
;   delete(zvarnew)
;   delete(wvarnew)
;   delete(utc_ztime)

; end reading Plev data
;*************************************************************************
; Ks data

;   fname = systemfunc("ls "+dirKs+"erain.Ks_wind.daily."+iyr+".nc")
;   ; fname = systemfunc("ls "+dirKs+"erain.Ks_wind.monmean."+iclyr+".nc")
;   print ("reading Ks file "+ fname)
;   Ksfile = addfile(fname,"r")

;   Kslev = Ksfile->level
;   iKslev = ind(doubletofloat(Kslev).eq.level)

;   if (getfilevartypes(Ksfile,"Ks").eq."short")then
;    Ksvar = short2flt(Ksfile->Ks(:,iKslev,:,:))
;   else
;    Ksvar = Ksfile->Ks(:,iKslev,:,:)
;   end if


;   lon = Ksfile->lon
;   lat = Ksfile->lat
;   Kslon = dimsizes(Ksvar&lon)
;   Kslat = dimsizes(Ksvar&lat)

; ;removing 29Feb in Plev
;   Kstime = Ksfile->time
;   utc_Kstime = cd_calendar(Kstime, -2)
  
;   t1 = ind(utc_Kstime.eq.ch28Feb)
;   t2 = ind(utc_Kstime.eq.ch1Mar)

;   Kstimenew = array_append_record(Kstime(:t1),Kstime(t2:),0)
;   Ksvarnew  = doubletofloat(array_append_record(Ksvar(:t1,:,:),Ksvar(t2:,:,:),0))
;   copy_VarMeta(Ksvar(:364,:,:),Ksvarnew)

;   delete(Kstime)
;   delete(Ksvar)
  
;   Kstime = Kstimenew
;   Ksvar = Ksvarnew

;   delete(Kstimenew)
;   delete(Ksvarnew)
;   delete(utc_Kstime)

; end reading Ks data
;*************************************************************************


  if (iyr.eq.syr)then
    t2marr      = new ((/nyr,365,nlatsrfc,nlonsrfc/),float)
    mslparr     = new ((/nyr,365,nlatsrfc,nlonsrfc/),float)
    srfctimearr = new ((/nyr,365/),typeof(srfctime))
;     zarr        = new ((/nyr,365,zlat,zlon/),float)
;     ztimearr    = new ((/nyr,365/),typeof(ztime))
;     warr        = new ((/nyr,365,zlat,zlon/),float)
;     Ksarr       = new ((/nyr,365,Kslat,Kslon/),float)
;     Kstimearr    = new ((/nyr,365/),typeof(Kstime))
  end if


  mslparr(iyr-syr,:,:,:) = mslpvar
  t2marr(iyr-syr,:,:,:)  = tvar
  srfctimearr(iyr-syr,:)  = srfctime
;   zarr(iyr-syr,:,:,:) = zvar
;   warr(iyr-syr,:,:,:)  = wvar
;   ztimearr(iyr-syr,:)  = ztime
;   Ksarr(iyr-syr,:,:,:)  = Ksvar
;   Kstimearr(iyr-syr,:)  = Kstime


  delete(tvar)
  delete(mslpvar)
  delete(srfctime)
;   delete(zvar)
;   delete(wvar)
;   delete(ztime)
;   delete(Ksvar)
;   delete(Kstime)


end do  ; year

; climatology
print("Climatology...")

iy1 = ind(yrs.eq.syr_cl)
iy2 = ind(yrs.eq.eyr_cl)

mslpclim = dim_avg_n_Wrap(mslparr(iy1:iy2,:,:,:), 0)
t2mclim = dim_avg_n_Wrap(t2marr(iy1:iy2,:,:,:), 0)
; zclim = dim_avg_n_Wrap(zarr(iy1:iy2,:,:,:), 0)
; wclim = dim_avg_n_Wrap(warr(iy1:iy2,:,:,:), 0)
; Ksclim = dim_avg_n_Wrap(Ksarr(iy1:iy2,:,:,:), 0)

;daily anomalies
print("Daily anomalies...")
amslp = mslparr
amslp = amslp@_FillValue
at2m = t2marr
at2m = at2m@_FillValue
; az = zvar
; az = az@_FillValue
; aw = warr
; aw = aw@_FillValue
; aKs = Ksarr
; aKs = aKs@_FillValue
do iyr = syr,eyr
    amslp(iyr-syr,:,:,:) = mslparr(iyr-syr,:,:,:) - mslpclim
    at2m(iyr-syr,:,:,:)   = t2marr(iyr-syr,:,:,:)  - t2mclim
    ; az(iyr-syr,:,:,:)    = zarr(iyr-syr,:,:,:)    - zclim
    ; aw(iyr-syr,:,:,:)    = warr(iyr-syr,:,:,:)    - wclim
    ; aKs(iyr-syr,:,:,:)   = Ksarr(iyr-syr,:,:,:)   - Ksclim
end do


; end
print("Data read")
;  END reading data

; reading T max/Tmin dates
print("Reading dates of extremes from files")
tmax  = new((/2,nreg,nmax/),float,-999)
tmin  = tmax

avmslp = new((/2,2,nreg,nmax,nlatsrfc,nlonsrfc/),typeof(amslp))
avamslp = avmslp
avt2m = new((/2,2,nreg,nmax,nlatsrfc,nlonsrfc/),typeof(at2m))
avat2m = avt2m
; avz = new((/2,2,nreg,nmax,zlat,zlon/),typeof(az))
; avaz = avz
; avw = new((/2,2,nreg,nmax,zlat,zlon/),typeof(aw))
; avaw = avw
; avKs = new((/2,2,nreg,nmax,Kslat,Kslon/),typeof(aKs))
; avaKs = avKs

do ireg =0, nreg-1
    fin0  = addfile(dirin+ "dates.Tmax."+ndav+"dclimanom."+"reg"+ireg+"."+ssn+syr0+"_"+eyr0+".nc"    , "r") 
    fin1  = addfile(dirin+ "dates.dtTmax."+ndav+"dclimanom."+"reg"+ireg+"."+ssn+syr0+"_"+eyr0+".nc"    , "r") 

    tmax(0,ireg,:) = fin0->timemax
    tmin(0,ireg,:) = fin0->timemin

    tmax(1,ireg,:) = fin1->timemax
    tmin(1,ireg,:) = fin1->timemin

    do is=0,1  ; sign of extreme: 0 <-> "+" / 1 <-> "-"
    do i=0,nmax -1

       do tr = 0,1
         print("reg="+reg(ireg)+"  sign="+sign(is)+"  "+"  nmax="+i+" trend = "+trend(tr))
        if(is.eq.0)then
          date=cd_calendar(tmax(tr,ireg,i) , -2)
          utcdate = cd_calendar(tmax(tr,ireg,i) , 0)        
        else if (is.eq.1)then
          date=cd_calendar(tmin(tr,ireg,i) , -2)
          utcdate = cd_calendar(tmin(tr,ireg,i) , 0)        
        end if
        end if

        print(utcdate)
        year =  tointeger(utcdate(0,0)) 
        ; year=1983
        ; utcdate = 19830108
        print(year)
        printVarSummary(srfctimearr(:,:))

        utcsrfc = cd_calendar(srfctimearr(year-syr,:),-2)
        tsrfc = ind(utcsrfc.eq.date)
        t1srfc = tsrfc - tadj1
        t2srfc = tsrfc - tadj2

        ; utcz = cd_calendar(ztimearr(year-syr,:),-2)
        ; tz = ind(utcz.eq.utcdate)
        ; t1z = tz - tadj1
        ; t2z = tz - tadj2
        
        ; utcKs = cd_calendar(Kstimearr(year-syr,:),-2)
        ; tKs = ind(utcKs.eq.utcdate)
        ; t1Ks = tKs - tadj1
        ; t2Ks = tKs - tadj2
        
        iyr = year-syr
        print(date)
        ; print(utcsrfc)
        print(tsrfc)
        print(t1srfc+"  "+t2srfc+"  "+iyr)

        if(t1srfc.ge.0.and.t2srfc.ge.0)then
          avmslp(is,tr,ireg,i,:,:)  = dim_avg_n(mslparr(iyr,t2srfc:t1srfc,:,:),0)
          avamslp(is,tr,ireg,i,:,:) = dim_avg_n(amslp(iyr,t2srfc:t1srfc,:,:),0)
          avt2m(is,tr,ireg,i,:,:)   = dim_avg_n(t2marr(iyr,t2srfc:t1srfc,:,:),0)
          avat2m(is,tr,ireg,i,:,:)  = dim_avg_n(at2m(iyr,t2srfc:t1srfc,:,:),0)
        else if(t1srfc.ge.0.and.t2srfc.lt.0)then
          t2srfc=365+t2srfc
          avmslp(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(mslparr(iyr-1,t2srfc:,:,:),mslparr(iyr,:t1,:,:),0),0)
          avamslp(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(amslp(iyr-1,t2srfc:,:,:),amslp(iyr,:t1,:,:),0),0)
          avt2m(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(t2marr(iyr-1,t2srfc:,:,:),t2marr(iyr,:t1,:,:),0),0)
          avat2m(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(at2m(iyr-1,t2srfc:,:,:),at2m(iyr,:t1,:,:),0),0)
        else if(t1srfc.lt.0.and.t2srfc.lt.0)then
          t1srfc=365+t1srfc
          t2srfc=365+t2srfc
          avmslp(is,tr,ireg,i,:,:)  = dim_avg_n(mslparr(iyr-1,t2srfc:t1srfc,:,:),0) 
          avamslp(is,tr,ireg,i,:,:) = dim_avg_n(amslp(iyr-1,t2srfc:t1srfc,:,:),0) 
          avt2m(is,tr,ireg,i,:,:)  = dim_avg_n(t2marr(iyr-1,t2srfc:t1srfc,:,:),0) 
          avat2m(is,tr,ireg,i,:,:) = dim_avg_n(at2m(iyr-1,t2srfc:t1srfc,:,:),0) 
        end if
        end if
        end if

        ; if(t1z.ge.0.and.t2z.ge.0)then
        ;   avz(is,tr,ireg,i,:,:)     = dim_avg_n(zarr(iyr,t2z:t1z,:,:),0)
        ;   avaz(is,tr,ireg,i,:,:)    = dim_avg_n(az(iyr,t2z:t1z,:,:),0)
        ;   avw(is,tr,ireg,i,:,:)     = dim_avg_n(warr(iyr,t2z:t1z,:,:),0)
        ;   avaw(is,tr,ireg,i,:,:)    = dim_avg_n(aw(iyr,t2z:t1z,:,:),0)
        ; else if(t1z.ge.0.and.t2z.lt.0)then
        ;   t2z=365+t2z
        ;   avz(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(zarr(iyr-1,t2z:,:,:),zarr(iyr,:t1z,:,:),0),0)
        ;   avaz(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(az(iyr-1,t2z:,:,:),az(iyr,:t1z,:,:),0),0)
        ;   avw(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(warr(iyr-1,t2z:,:,:),warr(iyr,:t1z,:,:),0),0)
        ;   avaw(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(aw(iyr-1,t2z:,:,:),aw(iyr,:t1z,:,:),0),0)
        ; else if(t1z.lt.0.and.t2z.lt.0)then
        ;   t1z=365+t1z
        ;   t2z=365+t2z
        ;   avz(is,tr,ireg,i,:,:)  = dim_avg_n(zarr(iyr-1,t2z:t1z,:,:),0) 
        ;   avaz(is,tr,ireg,i,:,:) = dim_avg_n(az(iyr-1,t2z:t1z,:,:),0) 
        ;   avw(is,tr,ireg,i,:,:)  = dim_avg_n(warr(iyr-1,t2z:t1z,:,:),0) 
        ;   avaw(is,tr,ireg,i,:,:) = dim_avg_n(aw(iyr-1,t2z:t1z,:,:),0) 
        ; end if
        ; end if
        ; end if

        ; if(t1Ks.ge.0.and.t2Ks.ge.0)then
        ;   avKs(is,tr,ireg,i,:,:)    = dim_avg_n(Ksarr(iyr,t2Ks:t1Ks,:,:),0)
        ;   avaKs(is,tr,ireg,i,:,:)   = dim_avg_n(aKs(iyr,t2Ks:t1Ks,:,:),0)
        ; else if(t1Ks.ge.0.and.t2Ks.lt.0)then
        ;   t2Ks=365+t2Ks
        ;   avKs(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(Ksarr(iyr-1,t2Ks:,:,:),Ksarr(iyr,:t1Ks,:,:),0),0)
        ;   avaKs(is,tr,ireg,i,:,:) = dim_avg_n(array_append_record(aKs(iyr-1,t2Ks:,:,:),aKs(iyr,:t1Ks,:,:),0),0)
        ; else if(t1Ks.lt.0.and.t2Ks.lt.0)then
        ;   t1Ks=365+t1Ks
        ;   t2Ks=365+t2Ks
        ;   avKs(is,tr,ireg,i,:,:)  = dim_avg_n(Ksarr(iyr-1,t2Ks:t1Ks,:,:),0) 
        ;   avaKs(is,tr,ireg,i,:,:) = dim_avg_n(aKs(iyr-1,t2Ks:t1Ks,:,:),0) 
        ; end if
        ; end if
        ; end if

       end do 
    end do
    end do

end do

;*************************************************************************
; Plotting
;*************************************************************************
plot = new(4,graphic)
do ireg =0, nreg-1
  do tr = 0,1  ;raw / detrended
    do is = 0,1  ; tmax/tmin
        ;surface composite
        fcmst_srfc = "Srfc_"+sign(is)+"."+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+trend(tr)+"."+ndav+"dav."+ssn+syr+"_"+eyr
        wks_srfc   = gsn_open_wks ("png",diro+fcmst_srfc)
        print(diro+fcmst_srfc)

        res = True

        res@gsnDraw  = False            
        res@gsnFrame = False
        res@gsnPolar   = "NH"
        res@mpGeophysicalLineThicknessF = 2.

        res@mpMinLatF          = 20    ; range to zoom in on
        res@mpMaxLatF          = 90
        res@cnFillOn           = True
        res@cnLinesOn          = False                   ; turn off contour lines
        res@cnLineLabelsOn     = False

        res@tiMainString      = ""
        res@gsnRightString     = "Treg: "+ndav+" daily average"

        sres                      = True              ; set up a second resource list
        sres@gsnDraw              = False             ; do not draw the plot
        sres@gsnFrame             = False             ; do not advance the frame
        sres@cnFillOn             = False
        ;    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
        ;    sres@cnLevels             = ispan(-5,35,5 )   ; set the contour levels

        cmpstmslp  = dim_avg_n(avmslp(is,tr,ireg,:,:,:),0)
        cmpstamslp = dim_avg_n(avamslp(is,tr,ireg,:,:,:),0)
        cmpstmslpstd  = dim_stddev_n(avmslp(is,tr,ireg,:,:,:),0)
        cmpstamslpstd  = dim_stddev_n(avamslp(is,tr,ireg,:,:,:),0)

        cmpstmslp!0 = "lat"
        cmpstmslp!1 = "lon"
        cmpstmslp&lat = latsrfc
        cmpstmslp&lon = lonsrfc
        copy_VarMeta(cmpstmslp, cmpstamslp)
        copy_VarMeta(cmpstmslp, cmpstmslpstd)
        copy_VarMeta(cmpstmslp, cmpstamslpstd)

        cmpstt2m  = dim_avg_n(avt2m(is,tr,ireg,:,:,:),0)
        cmpstat2m = dim_avg_n(avat2m(is,tr,ireg,:,:,:),0)
        cmpstt2mstd  = dim_stddev_n(avt2m(is,tr,ireg,:,:,:),0)
        cmpstat2mstd  = dim_stddev_n(avat2m(is,tr,ireg,:,:,:),0)

        cmpstt2m!0 = "lat"
        cmpstt2m!1 = "lon"
        cmpstt2m&lat = latsrfc
        cmpstt2m&lon = lonsrfc
        copy_VarMeta(cmpstt2m, cmpstat2m)
        copy_VarMeta(cmpstt2m, cmpstt2mstd)
        copy_VarMeta(cmpstt2m, cmpstat2mstd)


        res@gsnLeftString      = "Mslp"
        plot(0)  = gsn_csm_contour_map_polar (wks_srfc,cmpstmslp,res) 
        plot_ov  = gsn_csm_contour (wks_srfc,cmpstmslpstd,sres)    ; create the temperature plot 
        overlay(plot(0),plot_ov)

        res@gsnLeftString      = "T2m"
        plot(1)  = gsn_csm_contour_map_polar (wks_srfc,cmpstt2m,res) 
        plot_ov  = gsn_csm_contour (wks_srfc,cmpstt2mstd,sres)    ; create the temperature plot 
        overlay(plot(1),plot_ov)     

        res@gsnLeftString      = "Mslp daily anomalies"
        plot(2)  = gsn_csm_contour_map_polar (wks_srfc,cmpstamslp,res) 
        plot_ov  = gsn_csm_contour (wks_srfc,cmpstamslpstd,sres)    ; create the temperature plot 
        overlay(plot(2),plot_ov)

        res@gsnLeftString      = "T2m daily anomalies"
        plot(3)  = gsn_csm_contour_map_polar (wks_srfc,cmpstat2m,res) 
        plot_ov  = gsn_csm_contour (wks_srfc,cmpstat2mstd,sres)    ; create the temperature plot 
        overlay(plot(3),plot_ov)     


        ; Panel plot
        pres          = True       
        pres@gsnPanelMainString    = "Composites of Mslp and t2m for "+sign(is)+ " events in "+ \
                                       str_squeeze(reg(ireg))+ ", " + ssn   
        pres@gsnPanelMainFontHeightF = 0.01
        ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
        pres@gsnPanelTop     = 0.95
        gsn_panel(wks_srfc,plot,(/2,2/),pres)

        ; frame(wks_srfc)

    end do
  end do
end do

end