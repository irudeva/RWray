; 3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)
; See (38) of Takaya and Nakamura (2001, JAS)
;
; Used data:
;  Monthly-mean data of NCEP/NCAR reanalys 1
;   geopotential height (hgt: m) ; hgt.mon.mean.nc
;
;  Monthly climatology of NCEP/NCAR reanalys 1
;   geopotential height (hgt: m) ; hgt.mon.mean.nc
;   air tempeature (air: degC)   ; air.mon.mean.nc
;   zonal wind (uwnd: m/s)       ; uwnd.mon.mean.nc
;   meridional wind (vwnd: m/s)  ; vwnd.mon.mean.nc
;
; The data are available at
; http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;
; The period for calcuation is January, 1990 to December, 1990
; (you can modify this by changing fyear,fmon, lyear, lmon)
;
; The unit of level is [hPa]
;
; Basic state: monthly climatology
;
; Perturbation: deviation from climatology
;
; Calculating the flux for each month of input data.
;
; Note that
; (i) "p" (pressure) in Eq. (5.7) is dimensionless (pressure/1000mb)


; Output files are
; x-component:  TN2001-Fx.monthly.1990.nc
; y-component:  TN2001-Fy.monthly.1990.nc
; z-component:  TN2001-Fz.monthly.1990.nc
; QG stream function anomaly: psidev.monthly.1990.nc
; Brunt Vaisala frequency: NN.monthly.1990.nc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
; The first date of calculation
fyear = 1979
fmon = 1

; The last date of calculation
lyear = 2017
lmon = 12

;  direcotry name
diri = "/Users/irudeva/work/DATA/ERAint/Plev/"


; monthly climatology data
do iyr = fyear, lyear
  fname = systemfunc("ls "+diri+"erain.hgt_air_wind.monmean."+iyr+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"t").eq."short")then
   tvar = short2flt(clfile->t)
  else
   tvar = clfile->t
  end if
  tvar = tvar + 273.15
  ; tvar_yr(iyr-fyear,:,:,:) = tvar(level|:,latitude|:,longitude|:,time|:)

  if (getfilevartypes(clfile,"z").eq."short")then
   zvar = short2flt(clfile->z)
  else
   zvar = clfile->z
  end if
  ; zvar_ts = zvar(level|:,latitude|:,longitude|:,time|:)

  if (getfilevartypes(clfile,"u").eq."short")then
   uvar = short2flt(clfile->u)
  else
   uvar = clfile->u
  end if
  ; uvar_ts = uvar(level|:,latitude|:,longitude|:,time|:)

  if (getfilevartypes(clfile,"v").eq."short")then
   vvar = short2flt(clfile->v)
  else
   vvar = clfile->v
  end if
  ; if (iyr.eq.2017) then
  ;   print(vvar)
  ;   status_exit(1)
  ; end if
  ; vvar_ts = vvar(level|:,latitude|:,longitude|:,time|:)

  if (iyr.eq.fyear) then

   zvararr = zvar
   tvararr = tvar
   uvararr = uvar
   vvararr = vvar

  else
   zvararr1 = zvararr
   tvararr1 = tvararr
   uvararr1 = uvararr
   vvararr1 = vvararr

   delete(zvararr)
   delete(tvararr)
   delete(uvararr)
   delete(vvararr)

   zvararr = array_append_record (zvararr1,zvar, 0)
   tvararr = array_append_record (tvararr1,tvar, 0)
   uvararr = array_append_record (uvararr1,uvar, 0)
   vvararr = array_append_record (vvararr1,vvar, 0)

   delete(zvararr1)
   delete(tvararr1)
   delete(uvararr1)
   delete(vvararr1)
  end if


end do

;monthly climatology
zvararr_ts = zvararr(level|:,latitude|:,longitude|:,time|:)
tvararr_ts = tvararr(level|:,latitude|:,longitude|:,time|:)
uvararr_ts = uvararr(level|:,latitude|:,longitude|:,time|:)
vvararr_ts = vvararr(level|:,latitude|:,longitude|:,time|:)

z_clmon = clmMonLLLT( zvararr_ts )
t_clmon = clmMonLLLT( tvararr_ts )
u_clmon = clmMonLLLT( uvararr_ts )
v_clmon = clmMonLLLT( vvararr_ts )

;dates, lev and lat, lon
time = zvararr&time
level = zvararr&level
lat = zvararr&latitude
lon = zvararr&longitude

ntime = dimsizes(time)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlevel = dimsizes(level)

utc_date = cd_calendar(zvararr&time,0)
syear = tointeger(utc_date(0,0))
smon = tointeger(utc_date(0,1))

; new "time array" for output data
fstep = (fyear - syear)*12 + fmon - smon
lstep = (lyear - syear)*12 + lmon - smon
ctime = time(fstep:lstep)
nctime = dimsizes(ctime)

;selected arrays
szvar = new((/nctime,nlevel,nlat,nlon/),float)

clzvar = new((/nctime,nlevel,nlat,nlon/),float)
cltvar = new((/nctime,nlevel,nlat,nlon/),float)
cluvar = new((/nctime,nlevel,nlat,nlon/),float)
clvvar = new((/nctime,nlevel,nlat,nlon/),float)

do istep = 0,nctime-1
  iyear = tointeger(utc_date(fstep + istep,0))
  imon = tointeger(utc_date(fstep + istep,1))

  clzvar(istep,:,:,:) = z_clmon(:,:,:,imon-1)
  cltvar(istep,:,:,:) = t_clmon(:,:,:,imon-1)
  cluvar(istep,:,:,:) = u_clmon(:,:,:,imon-1)
  clvvar(istep,:,:,:) = v_clmon(:,:,:,imon-1)

end do

printVarSummary(zvararr)

;height anomaly
szvar = zvararr(fstep:lstep,:,:,:)
zavar = szvar - clzvar

;  Gas constant
gc=290
;  Gravitational acceleration
ga=9.80665

;  Radius of the earth
re=6378388

; scale height
sclhgt=8000.

; pi
pi = atan(1.0)*4.

; Coriolis parameter
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0 = "lat"
f&lat = lat
f@_FillValue = zvararr@_FillValue

; missing for 10S - 10N
do ilat = 0, nlat-1
 if (abs(lat(ilat) ).lt. 10. ) then
  f(ilat)= f@_FillValue
 end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

; 1-D -> 4-D
printVarSummary(zvararr)
leveltmp = conform_dims(dimsizes(zvararr),level,1)
coslattmp = conform_dims(dimsizes(zvararr),coslat,2)
ftmp = conform_dims(dimsizes(zvararr),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(cltvar*(1000./leveltmp)^0.286,-sclhgt*log(level/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = cltvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(cluvar^2 + clvvar^2)
cumag@_FillValue = cluvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
; psidev = zavar*ga /ftmp
psidev = zavar /ftmp  ; as zvar is alreeady Hgt*g

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,2)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(level/1000),False,0,1)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(level/1000),False,0,1)

;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(level/1000),False,0,1)

xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

; Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask

;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cluvar/coslattmp * xuterm + clvvar * xvterm), cluvar.lt.5,False)

;y-component
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cluvar*xvterm + coslattmp*clvvar*yvterm),  cluvar.lt.5,False)

;z-component
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(cluvar*zuterm + coslattmp*clvvar*zvterm), cluvar.lt.5, False)

; for output
Fx!0 = "time"
Fx&time = ctime
Fx!1 = "level"
Fx&level = level
Fx!2 = "lat"
Fx&lat = lat
Fx!3 = "lon"
Fx&lon = lon

Fy!0 = "time"
Fy&time = ctime
Fy!1 = "level"
Fy&level = level
Fy!2 = "lat"
Fy&lat = lat
Fy!3 = "lon"
Fy&lon = lon

Fz!0 = "time"
Fz&time = ctime
Fz!1 = "level"
Fz&level = level
Fz!2 = "lat"
Fz&lat = lat
Fz!3 = "lon"
Fz&lon = lon

psidev!0 = "time"
psidev&time = ctime
psidev!1 = "level"
psidev&level = level
psidev!2 = "lat"
psidev&lat = lat
psidev!3 = "lon"
psidev&lon = lon

zdev = zavar/ga
zdev!0 = "time"
zdev&time = ctime
zdev!1 = "level"
zdev&level = level
zdev!2 = "lat"
zdev&lat = lat
zdev!3 = "lon"
zdev&lon = lon

NN!0 = "time"
NN&time = ctime
NN!1 = "level"
NN&level = level
NN!2 = "lat"
NN&lat = lat
NN!3 = "lon"
NN&lon = lon

copy_VarMeta(zdev, clzvar)
copy_VarMeta(Fx, cltvar)
copy_VarMeta(Fx, cluvar)
copy_VarMeta(Fx, clvvar)


Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"

psidev@units = "m^2/s"

zdev@units = "m**2 s**-2"

NN@units = "1/s^2"

dirout = "../output/WAF/"
ncWAF =     addfile(dirout+"WAF.mon"+fyear+"-"+lyear+".nc","c")
; ncFx =     addfile(dirout+"Fx.mon.nc","c")
; ncFy =     addfile(dirout+"Fy.mon.nc","c")
; ncFz =     addfile(dirout+"Fz.mon.nc","c")
; ncZdev   = addfile(dirout+"Zdev.mon.nc","c")
; ncpsidev = addfile(dirout+"psidev.mon.nc","c")
; ncNN =     addfile(dirout+"NN.mon.nc","c")

ncWAF->Fx = Fx
ncWAF->Fy = Fy
ncWAF->Fz = Fz
; ncWAF->psidev = psidev
ncWAF->Zdev = zdev
; ncFx->Fx = Fx
; ncFy->Fy = Fy
; ncFz->Fz = Fz
; ncpsidev->psidev = psidev
; ncZdev->Zdev = zdev

; ncNN->NN = NN

; ncWAF->clzvar=clzvar
; ncWAF->ctvar=cltvar
ncWAF->cluvar=cluvar
ncWAF->clvvar=clvvar


end
