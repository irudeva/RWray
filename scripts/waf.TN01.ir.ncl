; 3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)
; See (38) of Takaya and Nakamura (2001, JAS)
;
; Used data:
;  Daily-mean data of NCEP/NCAR reanalys 1
;   geopotential height (hgt: m) ; hgt.1999.nc hgt.2000.nc hgt.2001.nc
;
;  Daily climatology of NCEP/NCAR reanalys 1 (365 days)
;   geopotential height (hgt: m) ; hgt.day.ltm.nc
;   air tempeature (air: degK)   ; air.day.ltm.nc
;   zonal wind (uwnd: m/s)       ; uwnd.day.ltm.nc
;   meridional wind (vwnd: m/s)  ; vwnd.day.ltm.nc
;
; The data are available at
; http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;
; The period for calcuation is 1st January, 2000 to 31th December, 2000
; (see variables of fyear,fmon, fday, lyear, lmon, lday)
;
; The unit of level is [hPa]
;
; Basic state: daily climatology
;
; Perturbation: deviation from climatology for each calendar day
;
; Calculating the flux for each day of input data.
;
; Note that
; (i)  10-day low-pass filter is applied
;       (Assume that phase speed (Cp) is almost zero)
;
; (ii) "p" (pressure) in Eq. (5.7) is dimensionless (pressure/1000mb)
;
; (iii) The year 2000 is a leap year (366 days)
;
; (iv)  Using height data of 1999 and 2001 is for time filtering.
;       The filter makes the first and last several days unavailable.

; Output files are
; x-component:  TN2001-Fx.daily.2000.nc
; y-component:  TN2001-Fy.daily.2000.nc
; z-component:  TN2001-Fz.daily.2000.nc
; daily QG stream function anomaly: psidev.daily.2000.nc
; Brunt Vaisala frequency: NN.daily.2000.nc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
chmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

; The first date of calculation
fyear = 2015
fmon = 7
fday = 1

; The last date of calculation
lyear = 2015
lmon = 7
lday = 31

;  direcotry name
diri = "/Users/irudeva/work/DATA/ERAint/Plev/"
dirci = diri

dirout = "../output/WAF/"

;
; daily mean data
; modify direcotry name
; filenames = systemfunc("ls "+diri+"hgt."+fyear+".mon.nc "+diri+"hgt.2000.nc "+diri+"hgt.2001.nc")
; filenames = systemfunc("ls "+diri+"erain.hgt."+fyear+".nc "+diri+"hgt.2000.nc "+diri+"hgt.2001.nc")
filenames = systemfunc("ls "+diri+"erain.hgt.daily."+fyear+".nc "+diri+"hgt.2000.nc "+diri+"hgt.2001.nc")
zfile = addfiles(filenames,"r")

ListSetType (zfile, "cat")

; daily climatology data
filename = systemfunc("ls "+dirci+"erain.hgt_air_wind.mon."+fyear+".nc")
btfile = addfile(filename,"r")
filename = systemfunc("ls "+dirci+"erain.hgt_air_wind.mon."+fyear+".nc")
bzfile = addfile(filename,"r")
filename = systemfunc("ls "+dirci+"erain.hgt_air_wind.mon."+fyear+".nc")
bufile = addfile(filename,"r")
filename = systemfunc("ls "+dirci+"erain.hgt_air_wind.mon."+fyear+".nc")
bvfile = addfile(filename,"r")

;Select levels and time

time = zfile[0]->time
lat = zfile[0]->lat
lon = zfile[0]->lon
level = zfile[0]->level

flev = ind(level.eq.200)
llev = ind(level.eq.1000)
slevel = level(flev:llev)
; print(flev+" "+llev)

ntime = dimsizes(time)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlevel = dimsizes(slevel)

; for climatology
print("Set the time var")
btime = btfile->time
nbtime = dimsizes(btime)

; Calendar for input file
; Assume that the calendar is "standard" (gregorian).
time@calendar = "standard"
option = 0
option@calendar = time@calendar
utc_date = cd_calendar( time,option)

; corresponding date of time(0) (not necesarily the same as the first date of the data)
syear = utc_date(0,0)
smon = utc_date(0,1)
sday = utc_date(0,2)

; Calendar for input file of daily climatology
; Assume that the calendar is "standard" (Gregorian).
btime@calendar = "standard"
option = 0
option@calendar = btime@calendar
butc_date = cd_calendar( btime(0),option)

; The first date of btime (not necesarily the same as the first date of the data)
; bsyear = tointeger(butc_date(0,0))
; bsmon = tointeger(butc_date(0,1))
; bsday = tointeger(butc_date(0,2))
; bshour = tointeger(butc_date(0,3))
; bsmin = tointeger(butc_date(0,4))
; bssec = tointeger(butc_date(0,5))

; time interval of the calendar
dtime = time(1)-time(0)
dbtime = btime(1)-btime(0)

fhh = 0
fmin = 0
fsec = 0

; in time() (double)
fdate = cd_inv_calendar(fyear,fmon,fday,fhh,fmin,fsec,time@units,option)
fstep = tointeger((fdate - time(0))/dtime)

lhh = 18
lmin = 0
lsec = 0
ldate = cd_inv_calendar(lyear,lmon,lday,lhh,lmin,lsec,time@units,option)
lstep = tointeger((ldate - time(0))/dtime)

; new "time array" for output data
ctime = time(fstep:lstep)
nctime = dimsizes(ctime)


; Read data
; Test a variable's type prior to reading the variable
print ("Reading daily data...")
if (getfilevartypes(zfile[0],"z").eq."short")then
  zvarorg = short2flt(zfile[:]->z(fstep:lstep,flev:llev,:,:))
else
  zvarorg = zfile[:]->z(fstep:lstep,flev:llev,:,:)
end if
;  Gravitational acceleration
; ga=9.80665
; zvarorg = zvarorg/ga

print ("Reading climatological data...")
if (getfilevartypes(btfile,"t").eq."short")then
  btvar = short2flt(btfile->t(:,flev:llev,:,:))
else
  btvar = btfile->t(:,flev:llev,:,:)
end if

if (getfilevartypes(bzfile,"z").eq."short")then
  bzvar = short2flt(bzfile->z(:,flev:llev,:,:))
else
  bzvar = bzfile->z(:,flev:llev,:,:)
end if
; bzvar = bzvar/ga

if (getfilevartypes(bufile,"u").eq."short")then
  buvar = short2flt(bufile->u(:,flev:llev,:,:))
else
  buvar = bufile->u(:,flev:llev,:,:)
end if

if (getfilevartypes(bvfile,"v").eq."short")then
  bvvar = short2flt(bvfile->v(:,flev:llev,:,:))
else
  bvvar = bvfile->v(:,flev:llev,:,:)
end if


;10-day low-pass fileter
print ("Filtering...")
wgt = filwgts_lanczos(21, 0,0.1, -999.,1)
zvar = wgt_runave_n_Wrap(zvarorg,wgt,0,0)

; make arrays for climatology that are suitable for calculation
; bzvar <- bzvarorg
; btvar <- btvarorg

print ("making arrays for clim data")
czvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
ctvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
cuvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)
cvvar = new((/nctime,nlevel,nlat,nlon/),float,zvarorg@_FillValue)

czvar@_FillValue = zvarorg@_FillValue
ctvar@_FillValue = zvarorg@_FillValue
cuvar@_FillValue = zvarorg@_FillValue
cvvar@_FillValue = zvarorg@_FillValue

print ("Start integration")

do istep = 0,nctime-1

  print ("istep="+istep+" out of "+nctime)

  iyear = tointeger(utc_date(fstep + istep,0))
  imon = tointeger(utc_date(fstep + istep,1))+1
  iday = tointeger(utc_date(fstep + istep,2))

  print ("     "+iday+"/"+imon+"/"+iyear)

  ;The climatological value of leap day is replaced by 28th February
  if (isleapyear(iyear) .and. imon .eq. 2 .and. iday .eq. 29) then
  ; 29feb -> 28feb
    iday = 28
  end if

  ; for monthly mean climatological data
  iday = 1

  idate = cd_inv_calendar(iyear,imon,iday,fhh,fmin,fsec,btime@units,option)
  jstep  = tointeger((idate-btime(0))/dbtime)
  czvar(istep,:,:,:) = bzvar(jstep,:,:,:)
  ctvar(istep,:,:,:) = btvar(jstep,:,:,:)
  cuvar(istep,:,:,:) = buvar(jstep,:,:,:)
  cvvar(istep,:,:,:) = bvvar(jstep,:,:,:)

  print("     "+"istep="+istep+ "    jstep="+jstep)
end do

;height anomaly
printVarSummary(zvar)
printVarSummary(czvar)

print(zvar(10,3,100,100))
print(czvar(10,3,100,100))

zavar = zvar - czvar
print(zavar(10,3,100,100))

; delete(czvar)
; delete(bzvar)
; delete(btvar)
; delete(buvar)
; delete(bvvar)

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
print ("Creating 4d var")
leveltmp = conform_dims(dimsizes(zvar),slevel,1)
coslattmp = conform_dims(dimsizes(zvar),coslat,2)
ftmp = conform_dims(dimsizes(zvar),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(ctvar*(1000./leveltmp)^0.286,-sclhgt*log(slevel/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = ctvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(cuvar^2 + cvvar^2)
cumag@_FillValue = cuvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
; psidev = zavar*ga /ftmp
psidev = zavar /ftmp

print(ftmp(10,3,100,100))
print(psidev(10,3,100,100))
status_exit(1)

print ("Taking derivatives...")

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
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(slevel/1000),False,0,1)
;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(slevel/1000),False,0,1)
;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(slevel/1000),False,0,1)

xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)

;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

; Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask
print ("Masking")

;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.5,False)

;y-component
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),  cuvar.lt.5,False)

;z-component
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(cuvar*zuterm + coslattmp*cvvar*zvterm), cuvar.lt.5, False)

print ("Writing the output...")

; for output
Fx!0 = "time"
Fx&time = ctime
Fx!1 = "level"
Fx&level = slevel
Fx!2 = "lat"
Fx&lat = lat
Fx!3 = "lon"
Fx&lon = lon

Fy!0 = "time"
Fy&time = ctime
Fy!1 = "level"
Fy&level = slevel
Fy!2 = "lat"
Fy&lat = lat
Fy!3 = "lon"
Fy&lon = lon

Fz!0 = "time"
Fz&time = ctime
Fz!1 = "level"
Fz&level = slevel
Fz!2 = "lat"
Fz&lat = lat
Fz!3 = "lon"
Fz&lon = lon

psidev!0 = "time"
psidev&time = ctime
psidev!1 = "level"
psidev&level = slevel
psidev!2 = "lat"
psidev&lat = lat
psidev!3 = "lon"
psidev&lon = lon

NN!0 = "time"
NN&time = ctime
NN!1 = "level"
NN&level = slevel
NN!2 = "lat"
NN&lat = lat
NN!3 = "lon"
NN&lon = lon

Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"

psidev@units = "m^2/s"

NN@units = "1/s^2"

chtime = fday+chmon(fmon)+fyear +"_"+lday+chmon(lmon)+lyear

fncFx =    dirout+"TN2001-Fx.daily."+chtime+".nc"
fncFy =    dirout+"TN2001-Fy.daily."+chtime+".nc"
fncFz =    dirout+"TN2001-Fz.daily."+chtime+".nc"
fncpsidev = dirout+"psidev.daily."+chtime+".nc"
fncNN =     dirout+"NN.daily."+chtime+".nc"

; remove any pre-existing file
system ("rm "+fncFx)
system ("rm "+fncFy)
system ("rm "+fncFz)
system ("rm "+fncpsidev)
system ("rm "+fncNN)

ncFx =     addfile(fncFx,"c")
ncFy =     addfile(fncFy,"c")
ncFz =     addfile(fncFz,"c")
ncpsidev = addfile(fncpsidev,"c")
ncNN =     addfile(fncNN,"c")

print("Writing Fx")
ncFx->Fx = Fx
print("Writing Fy")
ncFy->Fy = Fy
print("Writing Fz")
ncFz->Fz = Fz
print("Writing psidev")
ncpsidev->psidev = psidev
print("zavar to ncpsidev...")
copy_VarMeta(psidev,zavar)
ncpsidev->zavar = zavar
print("zvar to ncpsidev...")
ncpsidev->zvar = zvar
print("czvar to ncpsidev...")
;copy_VarMeta(psidev,czvar)
ncpsidev->czvar = czvar
print("zvarorg to ncpsidev...")
copy_VarMeta(psidev,zvarorg)
ncpsidev->zvarorg=zvarorg
print("ftmp to ncpsidev...")
copy_VarMeta(psidev,ftmp)
ncpsidev->ftmp=ftmp


print("Writing NN")
ncNN->NN = NN
end
