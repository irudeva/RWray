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

;var = "srfc"
var = "plev"
; var = 'Ks'

; reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
reg = (/"SSW" /)
nreg = dimsizes(reg)

ssn = "JFM"
if(ssn.eq."JFM")then
  t1ssn = 0
  t2ssn = 90
end if
if(ssn.eq."JAS")then
  t1ssn = 182
  t2ssn = 273
end if
syr0 = 1979
eyr0 = 2017

syr = 1979
eyr = 2017
yrs = ispan(syr,eyr,1)
nyr = eyr-syr+1

syr_cl = 1980 
eyr_cl = 2001
nyrs_cl = eyr_cl - syr_cl +1

; number of days for running average
ndav = 15

; number of extreme events
nmax = 20

ga=9.80665 ;  Gravitational acceleration

level = 300
;number of days for composits
ndcmpst = 10
tadj1 = 0
tadj2 = ndcmpst+tadj1-1

dirin   = "../output/SSW/" 
diro    = "../output/SSW/composites/"
dirPlev = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"
dirKs   = "/Users/irudeva/work/Projects/RWray/output/Ks/Ks_nc/"

; sign = (/"tSSW","tmin"/)
; sign = (/"warm","cold"/)
; trend = (/"raw","detrended"/)

;------------------------------------------------
; read data
;------------------------------------------------
print("Composties for "+var+" in " +syr +"-"+eyr+".")
print("  Climatology from "+syr_cl+"-"+eyr_cl)
print("    Number of days in composites "+ndcmpst+" from day "+tadj2+" to day " +tadj1)
sleep(5)

do iyr = syr,eyr
 
    ch28Feb = iyr+"0228"
    ch1Mar  = iyr+"0301"

  ;  Surface data
    if (var.eq."srfc")then

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
  end if ; var

;*************************************************************************
;  Plev data

  if (var.eq."plev")then
    fname = systemfunc("ls "+dirPlev+"erain.hgt_air_wind.daily."+iyr+".nc")
    print ("reading file "+ fname)
    clfile = addfile(fname,"r")

    zlev = clfile->level
    icllev = ind(zlev.eq.level)


    if (getfilevartypes(clfile,"z").eq."short")then
    zvar = short2flt(clfile->z(:,icllev,:,:))
    else
    zvar = clfile->z(:,icllev,:,:)
    end if
    zvar = zvar/ga
    lonz = zvar&lon
    latz = zvar&lat

    nlonz = dimsizes(lonz)
    nlatz = dimsizes(latz)


    print("getting wind components")
    if (getfilevartypes(clfile,"u").eq."short")then
    uvar = short2flt(clfile->u(:,icllev,:,:))
    else
    uvar = clfile->u(:,icllev,:,:)
    end if
    u2 = uvar*uvar

    if (getfilevartypes(clfile,"v").eq."short")then
    vvar = short2flt(clfile->v(:,icllev,:,:))
    else
    vvar = clfile->v(:,icllev,:,:)
    end if
    v2 = vvar*vvar
    print("calculating the wind speed...")
    wvar = sqrt(u2+v2)
    delete(uvar)
    delete(vvar)
    delete(u2)
    delete(v2)


    wvar!0 = "time"
    wvar&time = clfile->time
    wvar!1 = "latitude"
    wvar&latitude = clfile->lat
    wvar!2 = "longitude"
    wvar&longitude = clfile->lon
    print("end calculating the wind speed")

  ;removing 29Feb in Plev
    ztime = clfile->time
    utc_ztime = cd_calendar(ztime, -2)
    
    t1 = ind(utc_ztime.eq.ch28Feb)
    t2 = ind(utc_ztime.eq.ch1Mar)

    ztimenew = array_append_record(ztime(:t1),ztime(t2:),0)
    zvarnew = array_append_record(zvar(:t1,:,:),zvar(t2:,:,:),0)
    wvarnew = array_append_record(wvar(:t1,:,:),wvar(t2:,:,:),0)
    delete(ztime)
    delete(zvar)
    delete(wvar)

    ztime = ztimenew
    zvar = zvarnew
    wvar = wvarnew

    delete(ztimenew)
    delete(zvarnew)
    delete(wvarnew)
    delete(utc_ztime)

  end if ; var

; end reading Plev data
; *************************************************************************
; Ks data
  if (var.eq."Ks")then

    fname = systemfunc("ls "+dirKs+"erain.Ks_wind.daily."+iyr+".nc")
    ; fname = systemfunc("ls "+dirKs+"erain.Ks_wind.monmean."+iclyr+".nc")
    print ("reading Ks file "+ fname)
    Ksfile = addfile(fname,"r")

    Kslev = Ksfile->level
    iKslev = ind(doubletofloat(Kslev).eq.level)

    if (getfilevartypes(Ksfile,"Ks").eq."short")then
    Ksvar = short2flt(Ksfile->Ks(:,iKslev,:,:))
    else
    Ksvar = Ksfile->Ks(:,iKslev,:,:)
    end if


    lon = Ksfile->lon
    lat = Ksfile->lat
    lonKs = Ksvar&lon
    latKs = Ksvar&lat
    nlonKs = dimsizes(lonKs)
    nlatKs = dimsizes(latKs)

  ;removing 29Feb in Plev
    Kstime = Ksfile->time
    utc_Kstime = cd_calendar(Kstime, -2)
    
    t1 = ind(utc_Kstime.eq.ch28Feb)
    t2 = ind(utc_Kstime.eq.ch1Mar)

    Kstimenew = array_append_record(Kstime(:t1),Kstime(t2:),0)
    Ksvarnew  = doubletofloat(array_append_record(Ksvar(:t1,:,:),Ksvar(t2:,:,:),0))
    copy_VarMeta(Ksvar(:364,:,:),Ksvarnew)

    delete(Kstime)
    delete(Ksvar)
    
    Kstime = Kstimenew
    Ksvar = Ksvarnew

    delete(Kstimenew)
    delete(Ksvarnew)
    delete(utc_Kstime)

  end if ; var

; end reading Ks data
;*************************************************************************


  if (iyr.eq.syr)then
    if (var.eq."srfc")then
      t2marr      = new ((/nyr,365,nlatsrfc,nlonsrfc/),float)
      mslparr     = new ((/nyr,365,nlatsrfc,nlonsrfc/),float)
      srfctimearr = new ((/nyr,365/),typeof(srfctime))
    end if  ; var
    if (var.eq."plev")then
      zarr        = new ((/nyr,365,nlatz,nlonz/),float)
      ztimearr    = new ((/nyr,365/),typeof(ztime))
      warr        = new ((/nyr,365,nlatz,nlonz/),float)
    end if  ; var
    if (var.eq."Ks")then
      Ksarr       = new ((/nyr,365,nlatKs,nlonKs/),float)
      Kstimearr    = new ((/nyr,365/),typeof(Kstime))
    end if  ; var
  end if

  if (var.eq."srfc")then
    mslparr(iyr-syr,:,:,:) = mslpvar
    t2marr(iyr-syr,:,:,:)  = tvar
    srfctimearr(iyr-syr,:)  = srfctime
    delete(tvar)
    delete(mslpvar)
    delete(srfctime)
  end if  ; var
  if (var.eq."plev")then
    zarr(iyr-syr,:,:,:) = zvar
    warr(iyr-syr,:,:,:)  = wvar
    ztimearr(iyr-syr,:)  = ztime
    delete(zvar)
    delete(wvar)
    delete(ztime)
  end if  ; var
  if (var.eq."Ks")then
    Ksarr(iyr-syr,:,:,:)  = Ksvar
    Kstimearr(iyr-syr,:)  = Kstime
    delete(Ksvar)
    delete(Kstime)

    ;percentage of Ks not missing
    rKsarr = new((/nyr,365,nlatKs/),float)
    do iyr = syr,eyr
      do id =0,364
        do il = 0,nlatKs-1
          if(all(ismissing(Ksarr(iyr-syr,id,il,:))))then
            rKsarr(iyr-syr,id,il) = 0.
          else
            rKsarr(iyr-syr,id,il) = 100.*int2flt(num(.not.ismissing(Ksarr(iyr-syr,id,il,:))))/int2flt(dimsizes(Ksarr(iyr-syr,id,il,:)))
            ; print(num(.not.ismissing(Ksarr(iyr-syr,id,il,:)))+"  "+dimsizes(Ksarr(iyr-syr,id,il,:))+" rKs = "+rKsarr(iyr-syr,id,il))
          end if
        end do
      end do
    end do
  end if  ; var

end do  ; year

; climatology
print("Climatology...")

iy1 = ind(yrs.eq.syr_cl)
iy2 = ind(yrs.eq.eyr_cl)

if (var.eq."srfc")then
  mslpclim = dim_avg_n_Wrap(mslparr(iy1:iy2,:,:,:), 0)
  t2mclim = dim_avg_n_Wrap(t2marr(iy1:iy2,:,:,:), 0)
end if  ; var
if (var.eq."plev")then
  zclim = dim_avg_n_Wrap(zarr(iy1:iy2,:,:,:), 0)
  wclim = dim_avg_n_Wrap(warr(iy1:iy2,:,:,:), 0)
end if  ; var
if (var.eq."Ks")then
  Ksclim = dim_avg_n_Wrap(Ksarr(iy1:iy2,:,:,:), 0)
  rKsclim = dim_avg_n_Wrap(rKsarr(iy1:iy2,:,:), 0)

  rKsssn = dim_avg_n_Wrap(rKsarr(iy1:iy2,t1ssn:t2ssn,:), (/0,1/))
end if  ; var

;daily anomalies
print("Daily anomalies...")
if (var.eq."srfc")then
  amslp = mslparr
  amslp = amslp@_FillValue
  at2m = t2marr
  at2m = at2m@_FillValue
end if  ; var
if (var.eq."plev")then
  az = zarr
  az = az@_FillValue
  aw = warr
  aw = aw@_FillValue
end if  ; var
if (var.eq."Ks")then
  aKs = Ksarr
  aKs = aKs@_FillValue
  arKs = rKsarr
  arKs = arKs@_FillValue
end if  ; var
do iyr = syr,eyr
  if (var.eq."srfc")then
      amslp(iyr-syr,:,:,:) = mslparr(iyr-syr,:,:,:) - mslpclim
      at2m(iyr-syr,:,:,:)   = t2marr(iyr-syr,:,:,:)  - t2mclim
    end if  ; var
    if (var.eq."plev")then
      az(iyr-syr,:,:,:)    = zarr(iyr-syr,:,:,:)    - zclim
      aw(iyr-syr,:,:,:)    = warr(iyr-syr,:,:,:)    - wclim
    end if  ; var
    if (var.eq."Ks")then
      aKs(iyr-syr,:,:,:)   = Ksarr(iyr-syr,:,:,:)   - Ksclim
      arKs(iyr-syr,:,:)  = rKsarr(iyr-syr,:,:)   - rKsclim
    end if  ; var
end do


print("Data read")
;  END reading data

; reading T max/Tmin dates
print("Reading dates of extremes from files")
tSSW  = new((/nreg,nmax/),float,-999)
; tmin  = tSSW

if (var.eq."srfc")then
  avmslp = new((/nreg,nmax,nlatsrfc,nlonsrfc/),typeof(amslp))
  avamslp = avmslp
  avt2m = new((/nreg,nmax,nlatsrfc,nlonsrfc/),typeof(at2m))
  avat2m = avt2m
end if  ; var
if (var.eq."plev")then
  avz = new((/nreg,nmax,nlatz,nlonz/),typeof(az))
  avaz = avz
  avw = new((/nreg,nmax,nlatz,nlonz/),typeof(aw))
  avaw = avw
end if  ; var
if (var.eq."Ks")then
  avKs = new((/nreg,nmax,nlatKs,nlonKs/),typeof(aKs))
  avaKs = avKs
  avrKs = new((/nreg,nmax,nlatKs/),typeof(arKs))
  avarKs = avrKs
end if  ; var

do ireg =0, nreg-1
    fin0  = addfile(dirin+ "SSWdates.erain."+syr0+"_"+eyr0+".nc"    , "r") 

    tSSW(ireg,:) = fin0->timeSSW


    do i=0,nmax -1

        print("reg="+reg(ireg)+"  nmax="+i)
        date=cd_calendar(tSSW(ireg,i) , -2)
        utcdate = cd_calendar(tSSW(ireg,i) , 0)        
        year =  tointeger(utcdate(0,0)) 
        iyr = year-syr
        print(date)

        if (var.eq."srfc")then
          utcsrfc = cd_calendar(srfctimearr(iyr,:),-2)
          tsrfc = ind(utcsrfc.eq.date)
          t1srfc = tsrfc - tadj1
          t2srfc = tsrfc - tadj2
          ; print(utcsrfc)
          print(t1srfc)
        end if  ; var

        if (var.eq."plev")then
          utcz = cd_calendar(ztimearr(iyr,:),-2)
          tz = ind(utcz.eq.date)
          t1z = tz - tadj1
          t2z = tz - tadj2
        end if  ; var
        
        if (var.eq."Ks")then
          utcKs = cd_calendar(Kstimearr(iyr,:),-2)
          tKs = ind(utcKs.eq.date)
          t1Ks = tKs - tadj1
          t2Ks = tKs - tadj2
        end if  ; var

        if (var.eq."srfc")then
          print(t1srfc+"  "+t2srfc+"  "+iyr)
          if(t1srfc.ge.0.and.t2srfc.ge.0)then
            avmslp(ireg,i,:,:)  = dim_avg_n(mslparr(iyr,t2srfc:t1srfc,:,:),0)
            avamslp(ireg,i,:,:) = dim_avg_n(amslp(iyr,t2srfc:t1srfc,:,:),0)
            avt2m(ireg,i,:,:)   = dim_avg_n(t2marr(iyr,t2srfc:t1srfc,:,:),0)
            avat2m(ireg,i,:,:)  = dim_avg_n(at2m(iyr,t2srfc:t1srfc,:,:),0)
          else if(t1srfc.ge.0.and.t2srfc.lt.0)then
            t2srfc=365+t2srfc
            avmslp(ireg,i,:,:) = dim_avg_n(array_append_record(mslparr(iyr-1,t2srfc:,:,:),mslparr(iyr,:t1,:,:),0),0)
            avamslp(ireg,i,:,:) = dim_avg_n(array_append_record(amslp(iyr-1,t2srfc:,:,:),amslp(iyr,:t1,:,:),0),0)
            avt2m(ireg,i,:,:) = dim_avg_n(array_append_record(t2marr(iyr-1,t2srfc:,:,:),t2marr(iyr,:t1,:,:),0),0)
            avat2m(ireg,i,:,:) = dim_avg_n(array_append_record(at2m(iyr-1,t2srfc:,:,:),at2m(iyr,:t1,:,:),0),0)
          else if(t1srfc.lt.0.and.t2srfc.lt.0)then
            t1srfc=365+t1srfc
            t2srfc=365+t2srfc
            avmslp(ireg,i,:,:)  = dim_avg_n(mslparr(iyr-1,t2srfc:t1srfc,:,:),0) 
            avamslp(ireg,i,:,:) = dim_avg_n(amslp(iyr-1,t2srfc:t1srfc,:,:),0) 
            avt2m(ireg,i,:,:)  = dim_avg_n(t2marr(iyr-1,t2srfc:t1srfc,:,:),0) 
            avat2m(ireg,i,:,:) = dim_avg_n(at2m(iyr-1,t2srfc:t1srfc,:,:),0) 
          end if
          end if
          end if
        end if  ; var
        if (var.eq."plev")then
          print(t1z+"  "+t2z+"  "+iyr)
          if(t1z.ge.0.and.t2z.ge.0)then
            avz(ireg,i,:,:)     = dim_avg_n(zarr(iyr,t2z:t1z,:,:),0)
            avaz(ireg,i,:,:)    = dim_avg_n(az(iyr,t2z:t1z,:,:),0)
            avw(ireg,i,:,:)     = dim_avg_n(warr(iyr,t2z:t1z,:,:),0)
            avaw(ireg,i,:,:)    = dim_avg_n(aw(iyr,t2z:t1z,:,:),0)
          else if(t1z.ge.0.and.t2z.lt.0)then
            t2z=365+t2z
            avz(ireg,i,:,:) = dim_avg_n(array_append_record(zarr(iyr-1,t2z:,:,:),zarr(iyr,:t1z,:,:),0),0)
            avaz(ireg,i,:,:) = dim_avg_n(array_append_record(az(iyr-1,t2z:,:,:),az(iyr,:t1z,:,:),0),0)
            avw(ireg,i,:,:) = dim_avg_n(array_append_record(warr(iyr-1,t2z:,:,:),warr(iyr,:t1z,:,:),0),0)
            avaw(ireg,i,:,:) = dim_avg_n(array_append_record(aw(iyr-1,t2z:,:,:),aw(iyr,:t1z,:,:),0),0)
          else if(t1z.lt.0.and.t2z.lt.0)then
            t1z=365+t1z
            t2z=365+t2z
            avz(ireg,i,:,:)  = dim_avg_n(zarr(iyr-1,t2z:t1z,:,:),0) 
            avaz(ireg,i,:,:) = dim_avg_n(az(iyr-1,t2z:t1z,:,:),0) 
            avw(ireg,i,:,:)  = dim_avg_n(warr(iyr-1,t2z:t1z,:,:),0) 
            avaw(ireg,i,:,:) = dim_avg_n(aw(iyr-1,t2z:t1z,:,:),0) 
          end if
          end if
          end if
        end if  ; var

        if (var.eq."Ks")then
          print(t1Ks+"  "+t2Ks+"  "+iyr)
          if(t1Ks.ge.0.and.t2Ks.ge.0)then
            avKs(ireg,i,:,:)    = dim_avg_n(Ksarr(iyr,t2Ks:t1Ks,:,:),0)
            avaKs(ireg,i,:,:)   = dim_avg_n(aKs(iyr,t2Ks:t1Ks,:,:),0)
            avrKs(ireg,i,:)    = dim_avg_n(rKsarr(iyr,t2Ks:t1Ks,:),0)
            avarKs(ireg,i,:)   = dim_avg_n(arKs(iyr,t2Ks:t1Ks,:),0)
          else if(t1Ks.ge.0.and.t2Ks.lt.0)then
            t2Ks=365+t2Ks
            avKs(ireg,i,:,:) = dim_avg_n(array_append_record(Ksarr(iyr-1,t2Ks:,:,:),Ksarr(iyr,:t1Ks,:,:),0),0)
            avaKs(ireg,i,:,:) = dim_avg_n(array_append_record(aKs(iyr-1,t2Ks:,:,:),aKs(iyr,:t1Ks,:,:),0),0)
            avrKs(ireg,i,:) = dim_avg_n(array_append_record(rKsarr(iyr-1,t2Ks:,:),rKsarr(iyr,:t1Ks,:),0),0)
            avarKs(ireg,i,:) = dim_avg_n(array_append_record(arKs(iyr-1,t2Ks:,:),arKs(iyr,:t1Ks,:),0),0)
          else if(t1Ks.lt.0.and.t2Ks.lt.0)then
            t1Ks=365+t1Ks
            t2Ks=365+t2Ks
            avKs(ireg,i,:,:)  = dim_avg_n(Ksarr(iyr-1,t2Ks:t1Ks,:,:),0) 
            avaKs(ireg,i,:,:) = dim_avg_n(aKs(iyr-1,t2Ks:t1Ks,:,:),0) 
            avrKs(ireg,i,:)  = dim_avg_n(rKsarr(iyr-1,t2Ks:t1Ks,:),0) 
            avarKs(ireg,i,:) = dim_avg_n(arKs(iyr-1,t2Ks:t1Ks,:),0) 
          end if
          end if
          end if
        end if  ; var

    end do 


end do

;*************************************************************************
; Plotting
;*************************************************************************
plot = new(4,graphic)
do ireg =0, nreg-1
       ;surface composite
        if (var.eq."srfc")then
          fcmst_srfc = "Srfc_"+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+ssn+syr+"_"+eyr
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

          cmpstmslp  = dim_avg_n(avmslp(ireg,:,:,:),0)
          cmpstamslp = dim_avg_n(avamslp(ireg,:,:,:),0)
          cmpstmslpstd  = dim_stddev_n(avmslp(ireg,:,:,:),0)
          cmpstamslpstd  = dim_stddev_n(avamslp(ireg,:,:,:),0)

          cmpstmslp!0 = "lat"
          cmpstmslp!1 = "lon"
          cmpstmslp&lat = latsrfc
          cmpstmslp&lon = lonsrfc
          copy_VarMeta(cmpstmslp, cmpstamslp)
          copy_VarMeta(cmpstmslp, cmpstmslpstd)
          copy_VarMeta(cmpstmslp, cmpstamslpstd)

          cmpstt2m  = dim_avg_n(avt2m(ireg,:,:,:),0)
          cmpstat2m = dim_avg_n(avat2m(ireg,:,:,:),0)
          cmpstt2mstd  = dim_stddev_n(avt2m(ireg,:,:,:),0)
          cmpstat2mstd  = dim_stddev_n(avat2m(ireg,:,:,:),0)

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

          tres          = True
          title = "Composites of Mslp and t2m" 
          tres@txFont   = "helvetica"
          tres@txFontHeightF = 0.015
          gsn_text_ndc(wks_srfc,title,0.5,0.92,tres)

          ; Panel plot
          pres          = True       
          pres@gsnPanelMainString    = "Composites of Mslp and t2m "
          pres@gsnPanelMainFontHeightF = 0.01
          ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
          pres@gsnPanelTop     = 0.95
          gsn_panel(wks_srfc,plot,(/2,2/),pres)
        end if  ; var

        ; Plev compostes
        if (var.eq."plev")then
          fcmst_plev = "Plev_"+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+ssn+syr+"_"+eyr
          wks_plev   = gsn_open_wks ("png",diro+fcmst_plev)
          print(diro+fcmst_plev)

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

          cmpstz  = dim_avg_n(avz(ireg,:,:,:),0)
          cmpstaz = dim_avg_n(avaz(ireg,:,:,:),0)
          cmpstzstd  = dim_stddev_n(avz(ireg,:,:,:),0)
          cmpstazstd  = dim_stddev_n(avz(ireg,:,:,:),0)

          cmpstz!0 = "lat"
          cmpstz!1 = "lon"
          cmpstz&lat = latz
          cmpstz&lon = lonz
          copy_VarMeta(cmpstz, cmpstaz)
          copy_VarMeta(cmpstz, cmpstzstd)
          copy_VarMeta(cmpstz, cmpstazstd)

          cmpstw  = dim_avg_n(avw(ireg,:,:,:),0)
          cmpstaw = dim_avg_n(avaw(ireg,:,:,:),0)
          cmpstwstd  = dim_stddev_n(avw(ireg,:,:,:),0)
          cmpstawstd  = dim_stddev_n(avaw(ireg,:,:,:),0)

          cmpstw!0 = "lat"
          cmpstw!1 = "lon"
          cmpstw&lat = latz
          cmpstw&lon = lonz
          copy_VarMeta(cmpstw, cmpstaw)
          copy_VarMeta(cmpstw, cmpstwstd)
          copy_VarMeta(cmpstw, cmpstawstd)


          res@gsnLeftString      = "Hgt "+level
          plot(0)  = gsn_csm_contour_map_polar (wks_plev,cmpstz,res) 
          plot_ov  = gsn_csm_contour (wks_plev,cmpstzstd,sres)    ; create the temperature plot 
          overlay(plot(0),plot_ov)

          res@gsnLeftString      = "Wind at " + level +" hPa"
          plot(1)  = gsn_csm_contour_map_polar (wks_plev,cmpstw,res) 
          plot_ov  = gsn_csm_contour (wks_plev,cmpstwstd,sres)    ; create the temperature plot 
          overlay(plot(1),plot_ov)     

          res@gsnLeftString      = "Hgt "+level+ " daily anomalies"
          plot(2)  = gsn_csm_contour_map_polar (wks_plev,cmpstaz,res) 
          plot_ov  = gsn_csm_contour (wks_plev,cmpstazstd,sres)    ; create the temperature plot 
          overlay(plot(2),plot_ov)

          res@gsnLeftString      = "Wind daily anomalies"
          plot(3)  = gsn_csm_contour_map_polar (wks_plev,cmpstaw,res) 
          plot_ov  = gsn_csm_contour (wks_plev,cmpstawstd,sres)    ; create the temperature plot 
          overlay(plot(3),plot_ov)     

          tres          = True
          title = "Composites of Hgt "+level+" and wind" 
          tres@txFont   = "helvetica"
          tres@txFontHeightF = 0.015
          gsn_text_ndc(wks_plev,title,0.5,0.92,tres)

          ; Panel plot
          pres          = True       
          pres@gsnPanelMainString    = "Composites of Hgt "+level+" and wind"   
          pres@gsnPanelMainFontHeightF = 0.01
          ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
          pres@gsnPanelTop     = 0.95
          gsn_panel(wks_plev,plot,(/2,2/),pres)
        end if  ; var

        ; Ks plot
        if (var.eq."Ks")then
          fcmst_Ks = "Ks_"+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+ssn+syr+"_"+eyr
          wks_Ks   = gsn_open_wks ("png",diro+fcmst_Ks)
          print(diro+fcmst_Ks)

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

          cmpstKs  = dim_avg_n(avKs(ireg,:,:,:),0)
          cmpstaKs = dim_avg_n(avaKs(ireg,:,:,:),0)
          cmpstKsstd  = dim_stddev_n(avKs(ireg,:,:,:),0)
          cmpstaKsstd  = dim_stddev_n(avKs(ireg,:,:,:),0)

          cmpstKs!0 = "lat"
          cmpstKs!1 = "lon"
          cmpstKs&lat = latKs
          cmpstKs&lon = lonKs
          copy_VarMeta(cmpstKs, cmpstaKs)
          copy_VarMeta(cmpstKs, cmpstKsstd)
          copy_VarMeta(cmpstKs, cmpstaKsstd)

          res@gsnLeftString      = "Ks "+level
          res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
          res@cnMinLevelValF  = 0.			; set the minimum contour level
          res@cnMaxLevelValF  = 15.			; set the maximum contour level
          res@cnLevelSpacingF = 1.
          plot(0)  = gsn_csm_contour_map_polar (wks_Ks,cmpstKs,res) 
          plot_ov  = gsn_csm_contour (wks_Ks,cmpstKsstd,sres)    ; create the temperature plot 
          overlay(plot(0),plot_ov)

          cmpstrKs = dim_avg_n(avrKs(ireg,:,:),0)

          plotKs      = new((/2,dimsizes(cmpstKs&lat)/),float)
          plotKs(0,:) = rKsssn
          plotKs(1,:) = cmpstrKs

          kres                   = True
          kres@gsnDraw  = False            
          kres@gsnFrame = False
          kres@gsnLeftString     = "Ks profile at " + level +" hPa"
          kres@trYReverse        = True                       ; reverse Y-axis
          ; kres@xyDashPatterns    = 15                         ; choose dash patterns

          kres@xyLineThicknesses = (/  1.0,   2.0/)          ; make second line thicker
          kres@xyLineColors      = (/"blue","red"/)          ; change line color
          plot(1)  = gsn_csm_xy (wks_Ks,plotKs,cmpstKs&lat,kres) 
          delete(plotKs)
          
          res@gsnLeftString      = "Ks "+level+ " daily anomalies"
          res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
          res@cnMinLevelValF  = -7.			; set the minimum contour level
          res@cnMaxLevelValF  = 7.			; set the maximum contour level
          res@cnLevelSpacingF = 1.
          plot(2)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaKs,res) 
          plot_ov  = gsn_csm_contour (wks_Ks,cmpstaKsstd,sres)    ; create the temperature plot 
          overlay(plot(2),plot_ov)

          ; res@gsnLeftString      = "Wind daily anomalies"
          ; plot(3)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaw,res) 
          ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstawstd,sres)    ; create the temperature plot 
          ; overlay(plot(3),plot_ov)     

          tres          = True
          title = "Composites of Ks at "+level+" hPa " 
          tres@txFont   = "helvetica"
          tres@txFontHeightF = 0.015
          gsn_text_ndc(wks_Ks,title,0.5,0.92,tres)

          ; Panel plot
          pres = True       
          ; pres@gsnFrame = False
          pres@gsnPanelMainString    = "Composites of Ks at "+level+" hPa "  
          pres@gsnPanelMainFontHeightF = 0.015

          ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
          pres@gsnPanelTop     = 0.9
          pres@gsnPanelYWhiteSpacePercent = 5
          pres@gsnPanelXWhiteSpacePercent = 5
          print(pres)
          gsn_panel(wks_Ks,plot,(/2,2/),pres)
        end if  ; var

end do

end