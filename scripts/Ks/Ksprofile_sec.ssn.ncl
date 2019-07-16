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

chmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

; var = "srfc"
; var = "plev"
var = "Ks"

reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
nreg = dimsizes(reg)
chreg = stringtochar(reg)
chlen = dimsizes(chreg(0,:))-1


; reg = (/"C_Europe " /)
; nreg = dimsizes(reg)
lonregl = (/ 5, 30, 70, 100,  80, 100, 235, 260, 30, 160/)
lonregr = (/25, 50, 90, 120, 100, 120, 255, 280, 70, 200/)
lonregc = (lonregl+lonregr)/2.
; creating sectors
secwidth = 40 ; in deg lon
nsec = 360/secwidth + 1   ; 1 stands for zoal mean
lonsecl = new((/nreg,nsec/),float)
lonsecr = lonsecl
lonsecl(:,0) = 0.
lonsecr(:,0) = 360.
do ireg = 0, nreg-1
  do isec = 1,nsec-1
    lonsecl(ireg,isec) = lonregl(ireg)+(isec-1)*secwidth
    lonsecr(ireg,isec) = lonsecl(ireg,isec)+secwidth
  end do
end do

lonsecl = where(lonsecl.lt.0,lonsecl+360.,lonsecl)
lonsecl = where(lonsecl.ge.360.,lonsecl-360.,lonsecl)
lonsecr = where(lonsecr.le.0,lonsecr+360.,lonsecr)
lonsecr = where(lonsecr.gt.360.,lonsecr-360.,lonsecr)

; ssn = "JAS"

do issn = 0,1

if(issn.eq.0)then
  ssn = "JFM"
  t1ssn = 0
  t2ssn = 90
end if
if(issn.eq.1)then
  ssn = "JAS"
  t1ssn = 182  ;testing
  t2ssn = 273
end if
syr0 = 1980
eyr0 = 2017

syr = 1979
eyr = 2017  ; testing
yrs = ispan(syr,eyr,1)
nyr = eyr-syr+1

syr_cl = 1980 ;1981
eyr_cl = 2010  ; testing
nyrs_cl = eyr_cl - syr_cl +1

; number of days for running average
ndav = 15

; number of extreme events
nmax = 8

ga=9.80665 ;  Gravitational acceleration

level = 300
;number of days for composits
ndcmpst = 30
tadj1 = 0
tadj2 = ndcmpst+tadj1-1

dirin   = "../output/Treg/Treg.ssn/"
diro    = "../output/Treg/composites/ssn/"
dirPlev = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"
dirKs   = "/Users/irudeva/work/Projects/RWray/output/Ks/Ks_nc/"

; sign = (/"tmax","tmin"/)
sign = (/"warm","cold"/)
trend = (/"raw","detrended"/)

;------------------------------------------------
; read data
;------------------------------------------------
print("Composties for "+var+" in " + ssn+" "+syr +"-"+eyr+".")
print("  Climatology from "+syr_cl+"-"+eyr_cl)
sleep(5)


do iyr = syr,eyr

    ch28Feb = iyr+"0228"
    ch1Mar  = iyr+"0301"

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

    if (iyr.eq.syr)then
      rKsarr = new((/nyr,365,nlatKs/),float)
      rKsarrsec = new((/nyr,365,nreg,nsec,nlatKs/),float)
    end if

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

    if (var.eq."Ks")then
      Ksarr       = new ((/nyr,365,nlatKs,nlonKs/),float)
      Kstimearr    = new ((/nyr,365/),typeof(Kstime))
    end if  ; var
  end if


  if (var.eq."Ks")then
    Ksarr(iyr-syr,:,:,:)  = Ksvar
    Kstimearr(iyr-syr,:)  = Kstime
    delete(Ksvar)
    delete(Kstime)
  end if  ;var
  ;===================================================================
  ; calculating profiles of non-missing Ks
  ;===================================================================
    print("Calculating rKs for year = "+iyr)
    do id =t1ssn,t2ssn
      ; print("day="+id)
      print("iyr "+iyr+"  day "+id)
      do il = 0,nlatKs-1   ; testing
        ; print("day="+id+"  lat="+latKs(il))
        if(all(ismissing(Ksarr(iyr-syr,id,il,:))))then
          rKsarr(iyr-syr,id,il) = 0.
        else
          rKsarr(iyr-syr,id,il) = 100.*int2flt(num(.not.ismissing(Ksarr(iyr-syr,id,il,:))))/int2flt(dimsizes(Ksarr(iyr-syr,id,il,:)))
          ; print(num(.not.ismissing(Ksarr(iyr-syr,id,il,:)))+"  "+dimsizes(Ksarr(iyr-syr,id,il,:))+" rKs = "+rKsarr(iyr-syr,id,il))
        end if
        do ireg = 0,nreg-1
          do isec = 0, nsec-1
            ; print("ireg = "+ireg+"   isec="+isec)
            ; print(lonsecl(ireg,isec)+"  "+lonsecr(ireg,isec))
            ; print("iyr "+iyr+"  day "+id+"  lat: "+lat(il)+"    reg: "+reg(ireg)+"  sec: "+isec)
            if(lonsecl(ireg,isec).lt.lonsecr(ireg,isec))then
              Kstmp = Ksarr(iyr-syr,id,il,{lonsecl(ireg,isec):lonsecr(ireg,isec)-.01})
            else
              ; print("reverse")
              ; print(lonsecl(ireg,isec))
              ; print(lonsecr(ireg,isec))
              Kstmp = array_append_record(Ksarr(iyr-syr,id,il,{lonsecl(ireg,isec):}),Ksarr(iyr-syr,id,il,{:lonsecr(ireg,isec)-.01}),0)
            end if
            ; print("size "+dimsizes(Kstmp))
            if(all(ismissing(Kstmp)))then
              rKsarrsec(iyr-syr,id,ireg,isec,il) = 0.
            else
              rKsarrsec(iyr-syr,id,ireg,isec,il) = 100.*int2flt(num(.not.ismissing(Kstmp)))/int2flt(dimsizes(Kstmp))
              ; print(num(.not.ismissing(Ksarr(iyr-syr,id,il,:)))+"  "+dimsizes(Ksarr(iyr-syr,id,il,:))+" rKs = "+rKsarr(iyr-syr,id,il))
            end if
            ; print(reg(ireg)+" "+rKsarrsec(iyr-syr,id,ireg,isec,il))
            delete(Kstmp)
          end do  ; sector
        end do  ;reg


        end do
      end do


  ;===================================================================
  ; writing NETCDF
  ;===================================================================
  diro = "../output/Ks/Ksprofile/"
  filo = "Ksprof."+ssn+iyr+".nc"             ; Output file - problem with region name

  system("/bin/rm -f "  + diro + filo  )    ; remove if exists

  fout_rKs  = addfile(diro + filo    , "c")  ; open output file

  setfileoption(fout_rKs,"DefineMode",True)

  ;===================================================================
  ; create global attributes of the file
  ;===================================================================
  fAtt               = True            ; assign file attributes
  ; fAtt@title         = "NCL Efficient Approach to netCDF Creation"
  ; 	fAtt@source_file   =  "original-file.nc"
  fAtt@creation_date = systemfunc ("date")
  fAtt@Content   = "profiles of non-white Ks"
  ; fAtt@dataset = dataset
  ; fAtt@source = "Butler et al., 2017: A sudden stratospheric warming compendium."

  fileattdef( fout_rKs, fAtt )            ; copy file attributes

  ;===================================================================

  dimNames = (/ "nreg", "nsec", "nlat", "chlen"/)
  dimSizes = (/  nreg ,  nsec ,  nlatKs, chlen /)
  dimUnlim = (/ False, False, False, False/)
  filedimdef(fout_rKs,dimNames,dimSizes,dimUnlim)

  ;===================================================================

  filevardef(fout_rKs,     "lat" ,typeof(lat),                    "nlat" )
  print(typeof(reg))
  printVarSummary(reg)
  filevardef(fout_rKs,     "reg" ,typeof(chreg),      (/"nreg","chlen"/) )
  filevardef(fout_rKs, "reglonE" ,typeof(lonregr),            (/"nreg"/) )
  filevardef(fout_rKs, "reglonW" ,typeof(lonregl),            (/"nreg"/) )
  filevardef(fout_rKs, "seclonE" ,typeof(lonsecr),     (/"nreg","nsec"/) )
  filevardef(fout_rKs, "seclonW" ,typeof(lonsecl),     (/"nreg","nsec"/) )
  filevardef(fout_rKs,     "rKs" ,typeof(rKsarrsec),(/"nreg","nsec","nlat"/) )

  ;===================================================================
  ; Copy attributes associated with each variable to the file
  ; All attributes associated with each variable will be copied.
  ;====================================================================

  filevarattdef(fout_rKs,"lat",latKs)
  chreg@long_name = "name of region"
  filevarattdef(fout_rKs,"reg",chreg)
  lonregl@long_name = "longitude of west boundary of region"
  lonregl@units     = "degrees_east"
  filevarattdef(fout_rKs,"reglonW",lonregl)
  lonregl@units     = "degrees_east"
  lonregl@long_name = "longitude of east boundary of region"
  filevarattdef(fout_rKs,"reglonE",lonregl)
  lonsecl@units     = "degrees_east"
  lonsecl@long_name = "longitude of west boundary of sector"
  filevarattdef(fout_rKs,"seclonW",lonsecl)
  lonsecr@units     = "degrees_east"
  lonsecr@long_name = "longitude of east boundary of ssector"
  filevarattdef(fout_rKs,"seclonE",lonsecr)

  rKsout = rKsarrsec(0,0,:,:,:)
  rKsout@units     = "%"
  rKsout@long_name = "fraction of non-white Ks"
  filevarattdef(fout_rKs,"rKs",rKsout)
  delete(rKsout)


  ;===================================================================

  setfileoption(fout_rKs,"DefineMode",False)

  ;===================================================================

  ; fout_time->region   = stringtochar(str_squeeze(reg(ireg))  )
  fout_rKs->lat         = (/lat/)
  fout_rKs->reg         = (/chreg(:,:chlen-1)/)
  fout_rKs->reglonW     = (/lonregl/)
  fout_rKs->reglonE     = (/lonregr/)
  fout_rKs->seclonW     = (/lonsecl/)
  fout_rKs->seclonE     = (/lonsecr/)
  fout_rKs->rKs = (/dim_avg_n_Wrap(rKsarrsec(iyr-syr,t1ssn:t2ssn,:,:,:),(/0/))/)
; status_exit(1)
end do  ; year



; climatology
; print("Climatology...")
;
; iy1 = ind(yrs.eq.syr_cl)
; iy2 = ind(yrs.eq.eyr_cl)
;
; if (var.eq."Ks")then
;   ; Ksclim = dim_avg_n_Wrap(Ksarr(iy1:iy2,:,:,:), 0)
;   ; rKsclim = dim_avg_n_Wrap(rKsarr(iy1:iy2,:,:), 0)
;
;   rKsssnclim = dim_avg_n_Wrap(rKsarr(iy1:iy2,t1ssn:t2ssn,:), (/0,1/))
;   rKsssnsecclim = new((/nreg,nsec,dimsizes(rKsssnclim)/),typeof(rKsssnclim))
;   do ireg = 0,nreg-1
;     do isec = 0, nsec-1
;       rKsssnsecclim(ireg,isec,:) = dim_avg_n_Wrap(rKsarrsec(iy1:iy2,t1ssn:t2ssn,ireg,isec,:), (/0,1/))
;     end do
;   end do
;
;   ; rKsssn = dim_avg_n_Wrap(rKsarr(iy1:iy2,t1ssn:t2ssn,:), (/1/))
;   ; rKsssnsec = new((/nreg,nsec,dimsizes(rKsssn)/),typeof(rKsssn))
;   ; do ireg = 0,nreg-1
;   ;   do isec = 0, nsec-1
;   ;     rKsssnsec(:,ireg,isec,:) = dim_avg_n_Wrap(rKsarrsec(iy1:iy2,t1ssn:t2ssn,ireg,isec,:), (/1/))
;   ;   end do
;   ; end do
; end if  ; var

;daily anomalies
; print("Daily anomalies...")
; if (var.eq."srfc")then
;   amslp = mslparr
;   amslp = amslp@_FillValue
;   at2m = t2marr
;   at2m = at2m@_FillValue
; end if  ; var
; if (var.eq."plev")then
;   az = zarr
;   az = az@_FillValue
;   aw = warr
;   aw = aw@_FillValue
; end if  ; var
; if (var.eq."Ks")then
;   aKs = Ksarr
;   aKs = aKs@_FillValue
;   arKs = rKsarr
;   arKs = arKs@_FillValue
; end if  ; var
; do iyr = syr,eyr
;   if (var.eq."srfc")then
;       amslp(iyr-syr,:,:,:) = mslparr(iyr-syr,:,:,:) - mslpclim
;       at2m(iyr-syr,:,:,:)   = t2marr(iyr-syr,:,:,:)  - t2mclim
;     end if  ; var
;     if (var.eq."plev")then
;       az(iyr-syr,:,:,:)    = zarr(iyr-syr,:,:,:)    - zclim
;       aw(iyr-syr,:,:,:)    = warr(iyr-syr,:,:,:)    - wclim
;     end if  ; var
;     if (var.eq."Ks")then
;       aKs(iyr-syr,:,:,:)   = Ksarr(iyr-syr,:,:,:)   - Ksclim
;       arKs(iyr-syr,:,:)  = rKsarr(iyr-syr,:,:)   - rKsclim
;     end if  ; var
; end do


print("Finished reading data")
;  END reading data



; reading T max/Tmin dates
; print("Reading dates of extremes from files")
; yrmax  = new((/2,nreg,nmax/),float,-999)
; yrmin  = yrmax
;
;
;
; do ireg =0, nreg-1
;     fin  = addfile(dirin+ "t2m."+"reg"+ireg+"."+ssn+syr0+"_"+eyr0+".nc"   , "r")
;
;     yrmax(0,ireg,:) = fin->yr_t2m_max(:nmax-1)
;     yrmin(0,ireg,:) = fin->yr_t2m_min(:nmax-1)
;
;     yrmax(1,ireg,:) = fin->yr_t2mdt_max(:nmax-1)
;     yrmin(1,ireg,:) = fin->yr_t2mdt_min(:nmax-1)
;
;     print("years of raw t max: "+reg(ireg)+"  "+yrmax(0,ireg,:))
;
;     ; rKsssncmpst = dim_avg_n_Wrap(rKsarr(selind,t1ssn:t2ssn,:), (/0,1/))
;     if(ireg.eq.0)then
;       rKsssnseccmpst = new((/2,2,nreg,nsec,dimsizes(rKsssnclim)/),typeof(rKsssnclim))
;     end if
;
; ;for testing!!!
;     ; yrmax(0,ireg,:) = 1981
;     ; yrmin(0,ireg,:) = 1980
;     ;
;     ; yrmax(1,ireg,:) = 1980
;     ; yrmin(1,ireg,:) = 1982
;
;
;     do is=0,1  ; sign of extreme: 0 <-> "+" / 1 <-> "-"
;
;        do tr = 0,1
;          print("reg="+reg(ireg)+"  sign="+sign(is)+"  "+" trend = "+trend(tr))
;
;          selind = new(dimsizes(yrmax(is,ireg,:)),integer,0) ; will assing 1979 for all missing years
;          ; selind = new(dimsizes(yrmax(is,ireg,:)),integer)
;
;         do in = 0,  dimsizes(yrmax(is,ireg,:))-1
;           if(is.eq.0)then
;             selind(in) = ind(yrs.eq.yrmax(tr,ireg,in) )
;           else if (is.eq.1)then
;             selind(in) = ind(yrs.eq.yrmin(tr,ireg,in) )
;           end if
;           end if
;         end do
;
;         print("t ext years "+yrmax(is,ireg,:))
;         print("selected years "+yrs(selind))
;
;         if (var.eq."Ks")then
;
;          do isec = 0, nsec-1
;            rKsssnseccmpst(is, tr, ireg,isec,:) = dim_avg_n_Wrap(rKsarrsec(selind,t1ssn:t2ssn,ireg,isec,:), (/0,1/))
;            ; print(rKsarrsec(selind,1,ireg,isec,28:52))
;          end do
;
;         end if  ; var
;
;         ; print(reg(ireg)+" "+is+" "+tr+" rKsssnseccmpst "+rKsssnseccmpst(0,0,ireg,0,28:52))
;
;         delete(selind)
;
;        end do
;     end do
;
; end do
;
; ; status_exit(1)
;
; ;*************************************************************************
; ; Plotting
; ;*************************************************************************
; plot = new(4,graphic)
; plot_prof = new(nsec,graphic)
; do ireg =0, nreg-1
;   print("plotting1  "+reg(ireg)+"  "+rKsssnseccmpst(0,0,ireg,0,28:52))
;   do tr = 0,1  ;raw / detrended
;     do is = 0,1  ; tmax/tmin
;         ; Ks plot
;         if (var.eq."Ks")then
;           ; fcmst_Ks = "Ks_"+sign(is)+"."+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+trend(tr)+"."+ndav+"dav."+ssn+syr+"_"+eyr
;           ; wks_Ks   = gsn_open_wks ("png",diro+fcmst_Ks)
;           ; print(diro+fcmst_Ks)
;           ;
;           ; res = True
;           ;
;           ; res@gsnDraw  = False
;           ; res@gsnFrame = False
;           ; res@gsnPolar   = "NH"
;           ; res@mpGeophysicalLineThicknessF = 2.
;           ;
;           ; res@mpMinLatF          = 20    ; range to zoom in on
;           ; res@mpMaxLatF          = 90
;           ; res@cnFillOn           = True
;           ; res@cnLinesOn          = False                   ; turn off contour lines
;           ; res@cnLineLabelsOn     = False
;           ;
;           ; res@tiMainString      = ""
;           ; res@gsnRightString     = "Treg: "+ndav+" daily average"
;           ;
;           ; sres                      = True              ; set up a second resource list
;           ; sres@gsnDraw              = False             ; do not draw the plot
;           ; sres@gsnFrame             = False             ; do not advance the frame
;           ; sres@cnFillOn             = False
;           ; ;    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
;           ; ;    sres@cnLevels             = ispan(-5,35,5 )   ; set the contour levels
;           ;
;           ; cmpstKs  = dim_avg_n(avKs(is,tr,ireg,:,:,:),0)
;           ; cmpstaKs = dim_avg_n(avaKs(is,tr,ireg,:,:,:),0)
;           ; cmpstKsstd  = dim_stddev_n(avKs(is,tr,ireg,:,:,:),0)
;           ; cmpstaKsstd  = dim_stddev_n(avKs(is,tr,ireg,:,:,:),0)
;           ;
;           ; cmpstKs!0 = "lat"
;           ; cmpstKs!1 = "lon"
;           ; cmpstKs&lat = latKs
;           ; cmpstKs&lon = lonKs
;           ; copy_VarMeta(cmpstKs, cmpstaKs)
;           ; copy_VarMeta(cmpstKs, cmpstKsstd)
;           ; copy_VarMeta(cmpstKs, cmpstaKsstd)
;           ;
;           ; res@gsnLeftString      = "Ks "+level
;           ; res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;           ; res@cnMinLevelValF  = 0.			; set the minimum contour level
;           ; res@cnMaxLevelValF  = 15.			; set the maximum contour level
;           ; res@cnLevelSpacingF = 1.
;           ; plot(0)  = gsn_csm_contour_map_polar (wks_Ks,cmpstKs,res)
;           ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstKsstd,sres)    ; create the temperature plot
;           ; overlay(plot(0),plot_ov)
;
; ; plotting profiles
;             fprofile_Ks = "Ks_"+sign(is)+"."+str_squeeze(reg(ireg))+"."+trend(tr)+"."+ssn+syr+"_"+eyr
;
;             wks_type = "png"
;             wks_type@wkWidth = 2500
;             wks_type@wkHeight = 2500
;             wks_Ksprofile   = gsn_open_wks (wks_type,diro+fprofile_Ks)
;
;             print("writing to "+diro+fprofile_Ks)
;
;             plotKs      = new((/2,nlatKs/),float)
;             ; plotKs(0,:) = rKsssnsecclim(ireg,isec,:)
;             ; plotKs(1,:) = cmpstrKs
;             plotKs!1 = "lat"
;             plotKs&lat = latKs
;
;             kres                   = True
;             kres@gsnDraw  = False
;             kres@gsnFrame = False
;             ; kres@gsnLeftString     = "Ks profile at " + level +" hPa"
;             kres@trYReverse        = True                       ; reverse Y-axis
;             ; kres@xyDashPatterns    = 15                         ; choose dash patterns
;
;             kres@xyLineThicknesses = (/  3.,   4.0/)          ; make second line thicker
;             kres@xyLineColors      = (/"blue","red"/)          ; change line color
;             ; plot(1)  = gsn_csm_xy (wks_Ks,plotKs,cmpstKs&lat,kres)
;             ; delete(plotKs)
;
;             do isec = 0, nsec-1
;               kres@gsnLeftString     = "sector "+isec+" ("+lonsecl(ireg,isec)+"-"+lonsecr(ireg,isec)+"E)"
;               plotKs(0,:) = (/rKsssnsecclim(ireg,isec,:)/)
;               plotKs(1,:) = (/rKsssnseccmpst(is, tr, ireg,isec,:)/)
;               ; print("plotting  "+reg(ireg)+"  "+rKsssnseccmpst(is, tr,ireg,isec,28:52))
;               plot_prof(isec)  = gsn_csm_xy (wks_Ksprofile,plotKs,plotKs&lat,kres)
;             end do
;
;           ; res@gsnLeftString      = "Ks "+level+ " daily anomalies"
;           ; res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;           ; res@cnMinLevelValF  = -7.			; set the minimum contour level
;           ; res@cnMaxLevelValF  = 7.			; set the maximum contour level
;           ; res@cnLevelSpacingF = 1.
;           ; plot(2)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaKs,res)
;           ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstaKsstd,sres)    ; create the temperature plot
;           ; overlay(plot(2),plot_ov)
;
;           ; res@gsnLeftString      = "Wind daily anomalies"
;           ; plot(3)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaw,res)
;           ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstawstd,sres)    ; create the temperature plot
;           ; overlay(plot(3),plot_ov)
;
;           tres          = True
;           ; title = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
;           ;                               str_squeeze(reg(ireg))+ ", " + ssn
;           ; tres@txFont   = "helvetica"
;           ; tres@txFontHeightF = 0.015
;           ; gsn_text_ndc(wks_Ks,title,0.5,0.92,tres)
;           ;
;           ; title = "Profiles of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
;           ;                               str_squeeze(reg(ireg))+ ", " + ssn
;           ; gsn_text_ndc(wks_Ksprofile,title,0.5,0.92,tres)
;           ; ; Panel plot
;           ; pres = True
;           ; ; pres@gsnFrame = False
;           ; pres@gsnPanelMainString    = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
;           ;                               str_squeeze(reg(ireg))+ ", " + ssn
;           ; pres@gsnPanelMainFontHeightF = 0.015
;           ;
;           ; ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
;           ; pres@gsnPanelTop     = 0.9
;           ; pres@gsnPanelYWhiteSpacePercent = 5
;           ; pres@gsnPanelXWhiteSpacePercent = 5
;           ; print(pres)
;           ; gsn_panel(wks_Ks,plot,(/2,2/),pres)
;
;             title = "Profiles of Ks "+" at "+level+" hPa for "+" extreme "+sign(is)+ " seasons in "+ \
;                                           str_squeeze(reg(ireg))+ ", " + ssn
;             tres@txFont   = "helvetica"
;             tres@txFontHeightF = 0.015
;             gsn_text_ndc(wks_Ksprofile,title,0.5,0.92,tres)
;             ; LeftString = tointeger(utcdate(0,2))+" "+chmon(tointeger(utcdate(0,1))-1)+" "+tointeger(utcdate(0,0))
;             ; gsn_text_ndc(wks_Ksprofile,LeftString,0.1,0.88,tres)
;             ; RightString = "sector "+isec+" ("+lonsecl(ireg,isec)+"-"+lonsecr(ireg,isec)+"E)"
;             ; gsn_text_ndc(wks_Ksprofile,RightString,0.9,0.88,tres)
;
;             ; Panel plot
;             pres = True
;             ; pres@gsnFrame = False
;             ; pres@gsnPanelMainString    = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
;             ;                               str_squeeze(reg(ireg))+ ", " + ssn
;             ; pres@gsnPanelMainFontHeightF = 0.015
;
;             ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
;             pres@gsnPanelTop     = 0.9
;             pres@gsnPanelYWhiteSpacePercent = 5
;             pres@gsnPanelXWhiteSpacePercent = 5
;
;             ; if(is.eq.0)then
;             ;   date=cd_calendar(tmax(tr,ireg,imax) , -2)
;             ;   utcdate = cd_calendar(tmax(tr,ireg,imax) , 0)
;             ; else if (is.eq.1)then
;             ;   date=cd_calendar(tmin(tr,ireg,imax) , -2)
;             ;   utcdate = cd_calendar(tmin(tr,ireg,imax) , 0)
;             ; end if
;             ; end if
;
;             print(pres)
;             gsn_panel(wks_Ksprofile,plot_prof,(/5,2/),pres)
;
;           end if  ; var
;
;     end do
;   end do
; end do

end do ; ssn

end
