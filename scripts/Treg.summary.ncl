begin
;===================================================================
; summary plot
;===================================================================

diro = "../output/Treg/" 

syr = 1980
eyr = 2017
yrs = ispan(syr,eyr,1)
nyr = eyr-syr+1

ntmax = new((/2,nyr/),float,-999)
ntmax = 0
ntmin = ntmax


ssn ="JFM"

; number of days for running average
ndav = 15

; number of extreme events
nmax = 15

reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
nreg = dimsizes(reg)


tmax  = new((/2,nreg,nmax/),float,-999)
tmin  = tmax


do ireg =0, nreg-1
    fin0  = addfile(diro+ "dates.Tmax."+ndav+"dclimanom."+"reg"+ireg+"."+ssn+syr+"_"+eyr+".nc"    , "r") 
    fin1  = addfile(diro+ "dates.dtTmax."+ndav+"dclimanom."+"reg"+ireg+"."+ssn+syr+"_"+eyr+".nc"    , "r") 

    tmax(0,ireg,:) = fin0->timemax
    tmin(0,ireg,:) = fin0->timemin

    tmax(1,ireg,:) = fin1->timemax
    tmin(1,ireg,:) = fin1->timemin

    do i=0,nmax -1

       do tr = 0,1
        utc_date=cd_calendar(tmax(tr,ireg,i) , 0)
        year =  tointeger(utc_date(0,0)) 
        ntmax(tr,year - syr) = ntmax(tr,year - syr)+1
  
        utc_date=cd_calendar(tmin(tr,ireg,i) , 0)
        year =  tointeger(utc_date(0,0)) 
        ntmin(tr,year - syr) = ntmin(tr,year - syr)+1
       end do 
    end do

end do

;===================================================================
; plotting
;===================================================================	

wks   = gsn_open_wks ("png",diro+"T_summary."+ndav+"dclimanom."+ssn+syr+"_"+eyr)      
print(diro+"T_summary."+ndav+"dclimanom."+ssn+syr+"_"+eyr)

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15
     sres@trXMinF = 0-1
     sres@trXMaxF = nyr
     sres@trYMinF = -nreg*1.7
     sres@trYMaxF = nreg*1.7
     sres@gsnDraw = True
     sres@gsnFrame = False
     sres@gsnYRefLine           = 0.              ; reference line   
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.3           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = ispan(0,nyr-1,1)
     sres@tmXBLabels        = yrs
     sres@tmXBLabelAngleF      = 90
     sres@tmXBLabelFontHeightF = 0.0105
     sres@tmXTLabelFontHeightF = 0.0105
     sres@tmYLLabelFontHeightF = 0.0125
     sres@tiMainFontHeightF = 0.019
     sres@tiMainFont = "helvetica"
     sres@tiMainString = "Cumulative number of extreme events in "+nreg+" regions of the NH " 
     sres@gsnStringFontHeightF = 0.015
     sres@gsnRightString = ssn 
     sres@gsnLeftString = ndav + "-day average" 
    ;  sres@tiYAxisString = "(~S~o~N~C)"
     sres@tiXAxisString = "time, years"
     sres@tiXAxisFontHeightF  = 0.015
        
    
     sres@gsnXYBarChartColors = (/"red"/)	
     plot1 = gsn_csm_xy(wks,fspan(-.2,nyr-1.2,nyr),ntmax(0,:),sres)	
     plot3 = gsn_csm_xy(wks,fspan(-.2,nyr-1.2,nyr),-ntmin(0,:),sres)	
     sres@gsnXYBarChartColors = (/"lightblue"/)					
     plot2 = gsn_csm_xy(wks,fspan(0.2,nyr-.8,nyr),ntmax(1,:),sres)	
     plot3 = gsn_csm_xy(wks,fspan(0.2,nyr-.8,nyr),-ntmin(1,:),sres)




     lbres                    = True          ; labelbar only resources
     lbres@vpWidthF           = 0.3           ; labelbar width
     lbres@vpHeightF          = 0.1           ; labelbar height
     lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
     lbres@lbFillColors       = (/"lightblue"/)
     lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
     lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     lbres@lbPerimOn          = False
     lbres@lgPerimColor 	 = "white"
     labels = (/"detrended"/)
     gsn_labelbar_ndc(wks,1,labels,0.52,0.23,lbres)	; draw right labelbar column
	
     lbres@lbFillColors       = (/"red"/)
     labels = (/"raw data "/)
     gsn_labelbar_ndc(wks,1,labels,0.17,0.23,lbres)	; draw left labelbar column


    ;************************************************
    ; calculate the regression coefficient (slope)    
    ;************************************************
    rcmax0    = regline(yrs-syr, ntmax(0,:))          ; slope           
    rcmax1    = regline(yrs-syr, ntmax(1,:))          ; slope           
    ; rc@units = "degK/day"         
    regrmax   = new ( (/2,nyr/), float, -99)
    x = yrs-syr
  
    regrmax(0,:) = rcmax0*x + rcmax0@yintercept                      
    regrmax(1,:) = rcmax1*x + rcmax1@yintercept           


   
   res                     = True                   ; plot mods desired
   res@gsnDraw = True
   res@gsnFrame = False
   res@vpWidthF = 0.7
   res@vpHeightF = 0.5
   res@vpXF = .15
   res@trXMinF = 0-1
   res@trXMaxF = nyr
   res@trYMinF = -nreg*1.7
   res@trYMaxF = nreg*1.7
   res@tmXBLabelsOn = False
   res@tmYLLabelsOn = False

   res@xyMarkLineModes     = (/"Lines","Lines"/)  ; choose which have markers
;    res@xyMarkers           = 16                     ; choose type of marker 
;    res@xyMarkerColor       = "red"                  ; Marker color
;    res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)
   res@xyDashPatterns      = (/0,0/)                      ; solid line 
   res@xyLineThicknesses   = (/4,4/)                ; set second line to 2
   res@xyLineColors        = (/"red","lightblue"/)           
   res@tmYLFormat          = "f"                    ; not necessary but nicer labels 

   res@tiMainString        = ""  ; title


   plot5  = gsn_csm_xy (wks,x,regrmax,res)           ; create plot

    rcmin0    = regline(yrs-syr, ntmin(0,:))          ; slope           
    rcmin1    = regline(yrs-syr, ntmin(1,:))          ; slope           
    ; rc@units = "degK/day"         
    regrmin   = new ( (/2,nyr/), float, -99)

    regrmin(0,:) = rcmin0*x + rcmin0@yintercept                      
    regrmin(1,:) = rcmin1*x + rcmin1@yintercept  
    plot6  = gsn_csm_xy (wks,x,-regrmin,res)           ; create plot


frame(wks)      
end