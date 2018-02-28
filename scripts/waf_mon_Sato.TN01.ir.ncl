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
fyear = 2012
fmon = 12

; The last date of calculation
lyear = 2012
lmon = 12

posyr = (/1984,2004,2005,2006,2007,2008,2009,2011/)
negyr = (/1980,1981,1983,1987,1988,1989,1995,1996,2003/)
climyr = ispan(1979,2010,1)


;  direcotry name
diri = "/Users/irudeva/work/DATA/ERAint/Plev/"
; dirci = "/e3b/ncep/plev/monthly-climate/nc/"
;
; monthly mean data
; modify direcotry name
filename = systemfunc("ls "+diri+"erain.hgt_w_T.mnth.Dec.nc")
zfile = addfile(filename,"r")

; monthly climatology data
filename = systemfunc("ls "+diri+"erain.hgt_w_T.mnth.Dec.nc")
clfile = addfile(filename,"r")
;filename = systemfunc("ls "+dirci+"hgt.mon.ltm.nc")
;bzfile = addfile(filename,"r")
;filename = systemfunc("ls "+dirci+"uwnd.mon.ltm.1968-1996.nc")
;bufile = addfile(filename,"r")
;filename = systemfunc("ls "+dirci+"vwnd.mon.ltm.1968-1996.nc")
;bvfile = addfile(filename,"r")

; Test a variable's type prior to reading the variable
if (getfilevartypes(zfile,"z").eq."short")then
zvarorg = short2flt(zfile->z)
else
zvarorg = zfile->z
end if

if (getfilevartypes(clfile,"t").eq."short")then
btvar = short2flt(clfile->t) + 273.15
else
btvar = clfile->t + 273.15
end if

if (getfilevartypes(clfile,"z").eq."short")then
bzvar = short2flt(clfile->z)
else
bzvar = clfile->z
end if

if (getfilevartypes(clfile,"u").eq."short")then
buvar = short2flt(clfile->u)
else
buvar = clfile->u
end if

if (getfilevartypes(clfile,"v").eq."short")then
bvvar = short2flt(clfile->v)
else
bvvar = clfile->v
end if

time = zfile->time
level = zfile->level
lat = zfile->latitude
lon = zfile->longitude

ntime = dimsizes(time)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlevel = dimsizes(level)

; for climatology
;btime = btfile->time
;nbtime = dimsizes(btime)

; Calendar for input file
; Assume that the calendar is "standard" (gregorian).
time@calendar = "standard"
option = 0
option@calendar = time@calendar
utc_date = cd_calendar( time,option)

; print (utc_date)
syear = tointeger(utc_date(0,0))
smon = tointeger(utc_date(0,1))

; Calendar for input file of climatology
; Assume that the calendar is "standard" (Gregorian).
;btime@calendar = "standard"
;option = 0
;option@calendar = btime@calendar

; 2015/12/22 added (Thanks to Dr. Lu Guoyang.)
;-------------------
;butc_date = cd_calendar( btime,option)
; first year and month of climatology data
;bsyear = tointeger(butc_date(0,0))
;bsmon = tointeger(butc_date(0,1))
;--------------------

;fstep = (fyear - syear)*12 + fmon - smon
;lstep = (lyear - syear)*12 + lmon - smon

; new "time array" for output data

;ctime = time(fstep:lstep)
;nctime = dimsizes(ctime)
;nctime = dimsizes(posyr)
ctime = time(2012-syear)
nctime = 1

;zvar = zvarorg(fstep:lstep,:,:,:)
zvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)

; make arrays for climatology that are suitable for calculation
; bzvar <- bzvarorg
; btvar <- btvarorg

czvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
ctvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
cuvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
cvvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)

czvar@_FillValue = zvarorg@_FillValue
ctvar@_FillValue = zvarorg@_FillValue
cuvar@_FillValue = zvarorg@_FillValue
cvvar@_FillValue = zvarorg@_FillValue

do istep = 0,nctime-1
zvar(istep,:,:,:)= dim_avg_n(zvarorg(posyr-syear,:,:,:),0)


czvar(istep,:,:,:) = dim_avg_n(bzvar(negyr-syear,:,:,:),0)
ctvar(istep,:,:,:) = dim_avg_n(btvar(climyr-syear,:,:,:),0)
cuvar(istep,:,:,:) = dim_avg_n(buvar(climyr-syear,:,:,:),0)
cvvar(istep,:,:,:) = dim_avg_n(bvvar(climyr-syear,:,:,:),0)
end do

;height anomaly
zavar = zvar - czvar
delete(czvar)
delete(bzvar)
delete(btvar)
delete(buvar)
delete(bvvar)

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
f@_FillValue = zvarorg@_FillValue

; missing for 10S - 10N
do ilat = 0, nlat-1
if (abs(lat(ilat) ).lt. 10. ) then
f(ilat)= f@_FillValue
end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

; 1-D -> 4-D
leveltmp = conform_dims(dimsizes(zvar),level,1)
coslattmp = conform_dims(dimsizes(zvar),coslat,2)
ftmp = conform_dims(dimsizes(zvar),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(ctvar*(1000./leveltmp)^0.286,-sclhgt*log(level/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = ctvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(cuvar^2 + cvvar^2)
cumag@_FillValue = cuvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
;psidev = zavar*ga /ftmp    ; if z is in [m]
psidev = zavar /ftmp        ; if z is in [m2*s-2]

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
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.5,False)

;y-component
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),  cuvar.lt.5,False)

;z-component
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(cuvar*zuterm + coslattmp*cvvar*zvterm), cuvar.lt.5, False)

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

NN!0 = "time"
NN&time = ctime
NN!1 = "level"
NN&level = level
NN!2 = "lat"
NN&lat = lat
NN!3 = "lon"
NN&lon = lon

zanom = zavar / ga
zanom!0 = "time"
zanom&time = ctime
zanom!1 = "level"
zanom&level = level
zanom!2 = "lat"
zanom&lat = lat
zanom!3 = "lon"
zanom&lon = lon



Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"

psidev@units = "m^2/s"
zanom@units = "m  "

NN@units = "1/s^2"

dirout = "../output/WAF/Sato/"
ncFx =     addfile(dirout+"Fx.nc","c")
ncFy =     addfile(dirout+"Fy.nc","c")
ncFz =     addfile(dirout+"Fz.nc","c")
ncpsidev = addfile(dirout+"psidev.nc","c")
ncNN =     addfile(dirout+"NN.nc","c")

ncFx->Fx = Fx
ncFy->Fy = Fy
ncFz->Fz = Fz
; ncpsidev->psidev = psidev
ncpsidev->psidev = zanom ; to get [m]

ncNN->NN = NN
end
