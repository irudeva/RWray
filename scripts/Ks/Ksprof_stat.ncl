begin

var = "Ks"
level = 300


;reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
;nreg = dimsizes(reg)

; number of extreme events
nmax = 8
dirin   = "../output/Treg/Treg.ssn/"
sign = (/"warm","cold"/)
trend = (/"raw","detrended"/)


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

  ; climatology:
  syr0 = 1980
  eyr0 = 2017 ; testing
  clyrs = ispan(syr0,eyr0,1)


  syr = 1979
  eyr = 2017  ; testing
  yrs = ispan(syr,eyr,1)
  nyr = eyr-syr+1


  diro = "../output/Ks/Ksprofile/"

  do iyr=syr, eyr
    print("reading "+iyr)
    filo = "Ksprof."+ssn+iyr+".nc"             ; Output file - problem with region name

    f_rKs  = addfile(diro + filo    , "r")  ; open output file

    reg = f_rKs->reg
    nreg = dimsizes(reg(:,0))

    lonregE = f_rKs->reglonE
    lonregW = f_rKs->reglonW

    lonsecE = f_rKs->seclonE
    lonsecW = f_rKs->seclonW

    if(iyr.eq.syr)then
      lat = f_rKs->lat
      tmp_rKs = f_rKs->rKs
      dim_rKs = dimsizes(tmp_rKs)
      nsec = dim_rKs(1)
      nlat = dim_rKs(2)
      rKs = new((/dim_rKs(0),dim_rKs(1),nyr,dim_rKs(2)/),typeof(tmp_rKs))
      rKs!2 = "year"
      delete(tmp_rKs)
      delete(dim_rKs)
    end if
    rKs(:,:,iyr - syr,:) = f_rKs->rKs
    printVarSummary(rKs)

  ; status_exit(1)
  end do   ; year
  print("data read")

  ;Climatology
  print("Climatology...")
  ; print(get1Dindex(yrs,clyrs))
  rKsclim = dim_avg_n_Wrap(rKs(:,:,get1Dindex(yrs,clyrs),:),2)


  ; reading T max/Tmin dates
  print("Reading dates of extremes from files")
  yrmax  = new((/2,nreg,nmax/),float,-999)
  yrmin  = yrmax



  do ireg =0, nreg-1
      fin  = addfile(dirin+ "t2m."+"reg"+ireg+"."+ssn+syr0+"_"+eyr0+".nc", "r")
      ; print(fin)
      ; print("read temp from " + fin)

      yrmax(0,ireg,:) = fin->yr_t2m_max(:nmax-1)
      yrmin(0,ireg,:) = fin->yr_t2m_min(:nmax-1)

      yrmax(1,ireg,:) = fin->yr_t2mdt_max(:nmax-1)
      yrmin(1,ireg,:) = fin->yr_t2mdt_min(:nmax-1)

      print("years of raw t max: "+reg(ireg,:)+"  "+yrmax(0,ireg,:))

      if(ireg.eq.0)then
        rKscmpst = new((/2,2,nreg,nsec,nlat/),typeof(rKsclim))
      end if

  ;for testing!!!
      ; yrmax(0,ireg,:) = 1981
      ; yrmin(0,ireg,:) = 1980
      ;
      ; yrmax(1,ireg,:) = 1980
      ; yrmin(1,ireg,:) = 1982


      do is=0,1  ; sign of extreme: 0 <-> "+" / 1 <-> "-"

         do tr = 0,0
           print("reg="+reg(ireg,:)+"  sign="+sign(is)+"  "+" trend = "+trend(tr))
          ;
          ;  selind = new(dimsizes(yrmax(is,ireg,:)),integer,0) ; will assing 1979 for all missing years
          ;  ; selind = new(dimsizes(yrmax(is,ireg,:)),integer)
          ;
          ; do in = 0,  dimsizes(yrmax(is,ireg,:))-1
            if(is.eq.0)then
              ; selind(in) = ind(yrs.eq.yrmax(tr,ireg,in) )
              selind = get1Dindex(yrs,yrmax(tr,ireg,:))
            else if (is.eq.1)then
              ; selind(in) = ind(yrs.eq.yrmin(tr,ireg,in) )
              selind = get1Dindex(yrs,yrmin(tr,ireg,:))
            end if
            end if
          ; end do


          print("t ext years "+yrmax(is,ireg,:))
          print("selected years "+yrs(selind))

          if (var.eq."Ks")then

           do isec = 0, nsec-1
             rKscmpst(is, tr, ireg,isec,:) = dim_avg_n_Wrap(rKs(ireg,isec,selind,:), 0)
             ; print(rKsarrsec(selind,1,ireg,isec,28:52))
           end do

          end if  ; var

          ; print(reg(ireg)+" "+is+" "+tr+" rKsssnseccmpst "+rKsssnseccmpst(0,0,ireg,0,28:52))

          delete(selind)

        end do  ; tr
      end do  ; sign

    end do  ; reg

  ; status_exit(1)

  ;*************************************************************************
  ; Plotting
  ;*************************************************************************
  plot = new(4,graphic)
  plot_prof = new(nsec,graphic)
  do ireg =0, nreg-1
    ; print("plotting1  "+reg(ireg,:)+"  "+rKscmpst(0,0,ireg,0,28:52))
    do tr = 0,0  ;raw / detrended
      do is = 0,1  ; tmax/tmin
          ; Ks plot
          if (var.eq."Ks")then
            ; fcmst_Ks = "Ks_"+sign(is)+"."+str_squeeze(reg(ireg))+"."+ndcmpst+"dcmpst."+trend(tr)+"."+ndav+"dav."+ssn+syr+"_"+eyr
            ; wks_Ks   = gsn_open_wks ("png",diro+fcmst_Ks)
            ; print(diro+fcmst_Ks)
            ;
            ; res = True
            ;
            ; res@gsnDraw  = False
            ; res@gsnFrame = False
            ; res@gsnPolar   = "NH"
            ; res@mpGeophysicalLineThicknessF = 2.
            ;
            ; res@mpMinLatF          = 20    ; range to zoom in on
            ; res@mpMaxLatF          = 90
            ; res@cnFillOn           = True
            ; res@cnLinesOn          = False                   ; turn off contour lines
            ; res@cnLineLabelsOn     = False
            ;
            ; res@tiMainString      = ""
            ; res@gsnRightString     = "Treg: "+ndav+" daily average"
            ;
            ; sres                      = True              ; set up a second resource list
            ; sres@gsnDraw              = False             ; do not draw the plot
            ; sres@gsnFrame             = False             ; do not advance the frame
            ; sres@cnFillOn             = False
            ; ;    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
            ; ;    sres@cnLevels             = ispan(-5,35,5 )   ; set the contour levels
            ;
            ; cmpstKs  = dim_avg_n(avKs(is,tr,ireg,:,:,:),0)
            ; cmpstaKs = dim_avg_n(avaKs(is,tr,ireg,:,:,:),0)
            ; cmpstKsstd  = dim_stddev_n(avKs(is,tr,ireg,:,:,:),0)
            ; cmpstaKsstd  = dim_stddev_n(avKs(is,tr,ireg,:,:,:),0)
            ;
            ; cmpstKs!0 = "lat"
            ; cmpstKs!1 = "lon"
            ; cmpstKs&lat = latKs
            ; cmpstKs&lon = lonKs
            ; copy_VarMeta(cmpstKs, cmpstaKs)
            ; copy_VarMeta(cmpstKs, cmpstKsstd)
            ; copy_VarMeta(cmpstKs, cmpstaKsstd)
            ;
            ; res@gsnLeftString      = "Ks "+level
            ; res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
            ; res@cnMinLevelValF  = 0.			; set the minimum contour level
            ; res@cnMaxLevelValF  = 15.			; set the maximum contour level
            ; res@cnLevelSpacingF = 1.
            ; plot(0)  = gsn_csm_contour_map_polar (wks_Ks,cmpstKs,res)
            ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstKsstd,sres)    ; create the temperature plot
            ; overlay(plot(0),plot_ov)

  ; plotting profiles
              fprofile_Ks = "Ks_"+sign(is)+"."+str_squeeze(reg(ireg,:))+"."+trend(tr)+"."+ssn+syr+"_"+eyr

              wks_type = "png"
              wks_type@wkWidth = 2500
              wks_type@wkHeight = 2500
              wks_Ksprofile   = gsn_open_wks (wks_type,diro+"plot/"+fprofile_Ks)

              print("writing to "+diro+fprofile_Ks)

              plotKs      = new((/2,nlat/),float)
              ; plotKs(0,:) = rKsssnsecclim(ireg,isec,:)
              ; plotKs(1,:) = cmpstrKs
              plotKs!1 = "lat"
              plotKs&lat = lat

              kres                   = True
              kres@gsnDraw  = False
              kres@gsnFrame = False
              ; kres@gsnLeftString     = "Ks profile at " + level +" hPa"
              kres@trYReverse        = True                       ; reverse Y-axis
              ; kres@xyDashPatterns    = 15                         ; choose dash patterns

              kres@xyLineThicknesses = (/  3.,   4.0/)          ; make second line thicker
              kres@xyLineColors      = (/"blue","red"/)          ; change line color
              ; plot(1)  = gsn_csm_xy (wks_Ks,plotKs,cmpstKs&lat,kres)
              ; delete(plotKs)

              do isec = 0, nsec-1
                kres@gsnLeftString     = "sector "+isec+" ("+lonsecW(ireg,isec)+"-"+lonsecE(ireg,isec)+"E)"
                plotKs(0,:) = (/rKsclim(ireg,isec,:)/)
                plotKs(1,:) = (/rKscmpst(is, tr, ireg,isec,:)/)
                ; print("plotting  "+reg(ireg)+"  "+rKsssnseccmpst(is, tr,ireg,isec,28:52))
                plot_prof(isec)  = gsn_csm_xy (wks_Ksprofile,plotKs,plotKs&lat,kres)
              end do

            ; res@gsnLeftString      = "Ks "+level+ " daily anomalies"
            ; res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
            ; res@cnMinLevelValF  = -7.			; set the minimum contour level
            ; res@cnMaxLevelValF  = 7.			; set the maximum contour level
            ; res@cnLevelSpacingF = 1.
            ; plot(2)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaKs,res)
            ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstaKsstd,sres)    ; create the temperature plot
            ; overlay(plot(2),plot_ov)

            ; res@gsnLeftString      = "Wind daily anomalies"
            ; plot(3)  = gsn_csm_contour_map_polar (wks_Ks,cmpstaw,res)
            ; plot_ov  = gsn_csm_contour (wks_Ks,cmpstawstd,sres)    ; create the temperature plot
            ; overlay(plot(3),plot_ov)

            tres          = True
            ; title = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
            ;                               str_squeeze(reg(ireg))+ ", " + ssn
            ; tres@txFont   = "helvetica"
            ; tres@txFontHeightF = 0.015
            ; gsn_text_ndc(wks_Ks,title,0.5,0.92,tres)
            ;
            ; title = "Profiles of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
            ;                               str_squeeze(reg(ireg))+ ", " + ssn
            ; gsn_text_ndc(wks_Ksprofile,title,0.5,0.92,tres)
            ; ; Panel plot
            ; pres = True
            ; ; pres@gsnFrame = False
            ; pres@gsnPanelMainString    = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
            ;                               str_squeeze(reg(ireg))+ ", " + ssn
            ; pres@gsnPanelMainFontHeightF = 0.015
            ;
            ; ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
            ; pres@gsnPanelTop     = 0.9
            ; pres@gsnPanelYWhiteSpacePercent = 5
            ; pres@gsnPanelXWhiteSpacePercent = 5
            ; print(pres)
            ; gsn_panel(wks_Ks,plot,(/2,2/),pres)

              title = "Profiles of Ks "+" at "+level+" hPa for "+" extreme "+sign(is)+ " seasons in "+ \
                                            str_squeeze(reg(ireg,:))+ ", " + ssn
              tres@txFont   = "helvetica"
              tres@txFontHeightF = 0.015
              gsn_text_ndc(wks_Ksprofile,title,0.5,0.92,tres)
              ; LeftString = tointeger(utcdate(0,2))+" "+chmon(tointeger(utcdate(0,1))-1)+" "+tointeger(utcdate(0,0))
              ; gsn_text_ndc(wks_Ksprofile,LeftString,0.1,0.88,tres)
              ; RightString = "sector "+isec+" ("+lonsecl(ireg,isec)+"-"+lonsecr(ireg,isec)+"E)"
              ; gsn_text_ndc(wks_Ksprofile,RightString,0.9,0.88,tres)

              ; Panel plot
              pres = True
              ; pres@gsnFrame = False
              ; pres@gsnPanelMainString    = "Composites of Ks "+" at "+level+" hPa for extreme "+sign(is)+ " events in "+ \
              ;                               str_squeeze(reg(ireg))+ ", " + ssn
              ; pres@gsnPanelMainFontHeightF = 0.015

              ; pres@gsnPanelFigureStrings = (/"A","B","C"/)
              pres@gsnPanelTop     = 0.9
              pres@gsnPanelYWhiteSpacePercent = 5
              pres@gsnPanelXWhiteSpacePercent = 5

              ; if(is.eq.0)then
              ;   date=cd_calendar(tmax(tr,ireg,imax) , -2)
              ;   utcdate = cd_calendar(tmax(tr,ireg,imax) , 0)
              ; else if (is.eq.1)then
              ;   date=cd_calendar(tmin(tr,ireg,imax) , -2)
              ;   utcdate = cd_calendar(tmin(tr,ireg,imax) , 0)
              ; end if
              ; end if

              print(pres)
              gsn_panel(wks_Ksprofile,plot_prof,(/5,2/),pres)

            end if  ; var

      end do
    end do
  end do



  status_exit(1)

end do ; ssn
end
