;these files are loaded by default
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

fwnd = "wnd300.mnth.erain.nc"

;print(!pwd)
fw    = addfile(fwnd,"r")

time  = fw->time
;level = fw->level
lon   = fw->longitude
lat   = fw->latitude

ntime  = dimsizes(time)
;nlev   = dimsizes(level)
nlat   = dimsizes(lat)
nlon   = dimsizes(lon)

;; Select level
;ilev  = ind(level.eq.300)
; print("ilev="+ilev)

;;Select date - Sep2014
date = cd_calendar(time, 0)                  
idate = ind(date(:,0).eq.2014.and.date(:,1).eq.9)
 print("idate="+idate)
status_exit(0)

;u  = fu->uwnd(idate,ilev,:,:)
;v  = fv->vwnd(idate,ilev,:,:)
u  = fw->u(:,:,:)
v  = fw->v(:,:,:)
 printVarSummary(u)

;Stream function
;input rightmost two dimensions must be nlat x nlon!!!

 sfvp  = uv2sfvpF (u,v)
 printVarSummary(sfvp)
 sf    = sfvp(0,:,:,:)
  sf!0  = u!0
  sf!1  = u!1
  sf!2  = u!2
  sf@units = "m*m/s"
  sf@long_name = "Stream function"
  ;sf@time  = u@time
  sf@level = u@level
 printVarSummary(sf)

;Writing sf to netcdf

;------------------------------------------------------------

 ;fon = "sf"+level(ilev)+".nc"
 fon = "sf300.mnth.erain.nc"
 system("/bin/rm -f " + fon) ; remove if exists
 fo = addfile(fon, "c")

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
 setfileoption(fo,"DefineMode",True)

; create global attributes of the file
;===================================================================
 fAtt               = True            ; assign file attributes
 fAtt@title         = "NCL generated netCDF file"
 fAtt@source_file   = fwnd
;fAtt@creation_date = systemfunc ("date")
 fileattdef(fo, fAtt)      ; copy file attributes


 dimNames = (/"time",  getvardims(lat), getvardims(lon)/)
 dimSizes = (/ ntime,  nlat,  nlon /)
 dimUnlim = (/ False, False, False/)
 filedimdef(fo, dimNames, dimSizes, dimUnlim)

;===================================================================
; defining variables
;===================================================================


 filevardef(fo, "time", typeof(time), getvardims(time))
 filevarattdef(fo,"time", time)                   ; copy time attributes
 fo->time   = (/time/)

 ;filevardef(fo, "level", typeof(level), getvardims(level) )
 ;filevarattdef(fo,"level", level)                     ; copy lev attributes
 ;fo->lev    = (/level/)

 filevardef(fo, "lat", typeof(lat), getvardims(lat))
 filevarattdef(fo,"lat", lat)                     ; copy lat attributes
 fo->lat    = (/lat/)

 filevardef(fo, "lon", typeof(lon), getvardims(lon))
 filevarattdef(fo,"lon", lon)                     ; copy lon attributes
 fo->lon    = (/lon/)

 filevardef(fo, "sf", typeof(sf), getvardims(sf))
 filevarattdef(fo,"sf", sf)                   
 fo->sf      = (/sf/)

;status_exit(1)

end
