begin

diro = "../output/SIE/" 

syr = 1980
eyr = 2017
nyr = eyr-syr+1

ssn = "JFM"
if (ssn.eq."JFM") then
   mStrt = 1 ; -month start
   mLast = 3; -month last
   maxndays = 90
   nmon = 3
   mon = (/1,2,3/)
   mday = (/31,28,31/)
else if  (ssn.eq."JAS") then
   mStrt = 7 ; -month start
   mLast = 9 ; -month last
   maxndays = 92
   nmon = 3
   mon = (/7,8,9/)
   mday = (/31,31,30/)
end if
end if

;create an empty SIE array
units  = "hours since 1900-01-01 00:00:00"  
time=new(nyr*maxndays,integer)
time@units = units
time@calendar = "standard"
aSIE =new(nyr*maxndays,float)

h =0
m = 0
sec =0

i = 0
do iy = syr,eyr
   do im =  0,nmon-1
      do d = 1,mday(im)
        it = cd_inv_calendar(iy,mon(im),d,h,m,sec,units, 0)
        time(i) = tointeger(it)
        i = i+1
      end do
  end do
end do

; timeall=new((2018-1977+1)*364,integer)
; timeall@units = units
; timeall@calendar = "standard"
; aSIEall =new((2018-1977+1)*364,float)
; i = 0
; do iy = 1977,2018
;    do im =  0,11
;    nd = days_in_month(iy, im+1)
;    if(imon.eq.1)nd =28
;       do d = 1,nd
;         it = cd_inv_calendar(iy,mon(im),d,h,m,sec,units, 0)
;         timeall(i) = tointeger(it)
;         i = i+1
;       end do
;   end do
; end do


; reading data

   ncol = 5
   data = readAsciiTable("../data/N_seaice_extent_daily_v3.0.csv", ncol, "float", (/2,0/))

   rawSIE = data(:,3)
   timeSIE = tointeger(data(:,:2))
   delete(data)

   ; climatological day of year
   ncol = 8 
   data = readAsciiTable("../data/N_seaice_extent_climatology_1981-2010_v3.0.csv", ncol, "float", (/2,0/))

;    print(data(:,1))
   clSIE = data(:,1)
   clDOY = data(:,0)
   delete(data)   

;anomalous SIE

doy = day_of_year(timeSIE(:,0),timeSIE(:,1),timeSIE(:,2))
; print(doy)
arawSIE = rawSIE - clSIE(doy-1)


;raw time
   hh    = timeSIE(:,2)                     ; create arrays [*] of required size
   mn    = hh
   sc    = hh

   hh    = 0                                  ; array syntax
   mn    = 0
   sc    = 0
                                
   tmp   = cd_inv_calendar(timeSIE(:,0),timeSIE(:,1),timeSIE(:,2),hh,mn,sc,units, 0)
   rawtimeSIE = toint(tmp)
   rawtimeSIE!0 = "time"
   rawtimeSIE@units = units
   
   delete(tmp)
   delete(hh)
   delete(mn)
   delete(sc)

;subsetting

; seldate = ind(timeSIE(:,1).ge.mStrt.and.timeSIE(:,1).le.mLast.and.timeSIE(:,0).ge.syr.and.timeSIE(:,0).le.eyr)
; sel_aSIE = arawSIE(seldate)
; sel_hourSIE = hourSIE(seldate)

; writing to SIE
aSIE@_FillValue = arawSIE@_FillValue
do i=0,dimsizes(rawtimeSIE)-1
  if(.not.ismissing(ind(time.eq.rawtimeSIE(i))))then
   aSIE(ind(time.eq.rawtimeSIE(i))) = arawSIE(i)
  end if
end do

; print(aSIE(:100))
; status_exit(0)


;plotting

 fout = diro+"aSIE"+ssn+syr+"_"+eyr
  wks          = gsn_open_wks ("png",fout)

  res          = True                   ; plot mods desired
  res@gsnFrame = False                  ; don't advance frame yet      

  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.8                  
  
;   res@trYMinF  = -15.0                   ; min value on y-axis
;   res@trYMaxF  =  15.0                   ; max value on y-axis

  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012

  res@tmXBLabelAngleF      = 90

; since we stretch the plot, we need to slide it over a bit so it does not
; run off the page. we do this by:

  res@vpXF     = 0.1                   ; start plot at x ndc coord 

;   res@xyMarkLineMode = "Markers"
;   res@xyMarkers         =  16                      ; choose type of marker  
;   ;  res@xyMarkerColor     = "red"                    ; Marker color
;   res@xyMarkerSizeF     = 0.002                     ; Marker size (default 0.01)
  
  res@tiMainString     ="SIE, anomalies from climatological daily means, " + ssn 
  res@tiMainFontHeightF   = 0.015
;   res@gsnLeftString     = "window = "+window+" days"  
;   res@gsnStringFontHeightF   = 0.01
;   res@tiYAxisString    = "T anomalies, ~S~o~N~C"    ; y-axis label      
;   res@tiYAxisFontHeightF    = 0.012
;   res@tiXAxisString    = "time"    ; y-axis label      
;   res@tiXAxisFontHeightF    = 0.012
  
  y = aSIE
  x = ispan(1,dimsizes(aSIE),1)
;   x = time
;   y!0 = "time"
;   y&time = x

; interpolating y
y1 = linmsg_n(y, (/0,1/),0)


  yy = ispan(syr,eyr,1) 
  mm = yy
  dd = yy
  hh = yy
  mn = yy
  sc = yy

   mm = mon(0)
   dd = 1
   hh = 0                                 
   mn = 0
   sc = 0

   FirstJul = cd_inv_calendar(yy,mm,dd,hh,mn,sc,units, 0)

   Xval = new(nyr,integer)
   do i =0,nyr-1
    Xval(i)=ind(time.eq.FirstJul(i))
   end do

  res@tmXBMode        = "Explicit"	
  res@tmXBValues      = Xval  ;fspan(0,dimsizes(y),nyr+1)
  res@tmXBLabels      = "" + ispan(syr,eyr+1,1) 
  
  res@trXMinF  = min(x)-2                 
  res@trXMaxF  =  max(x)+2                

  plot = gsn_csm_xy (wks,x,y1,res)        ; create plot
;   gsn_polyline(wks,plot,dateF,(/dsoid/),polyres) ; add polyline
  frame(wks)                                     ; now advance frame

; print("x="+x)
; print("y="+y)
; printVarSummary(y)
; print(Xval)







end