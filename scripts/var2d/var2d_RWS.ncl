begin
;************************************************
; read in netCDF file
;************************************************

 ssn = "JJA"
 dir = "/users/Irina/work/Conference/AMOS2017/Poster/figs/"
 nvar = (/"u","v","um","vm","umx","umy","vmx","vmy","q","qx","qy","qxx","qyy","qxy","BetaM"/)
; do iv=0,dimsizes(nvar)-1
; do iv=0,1
  ;a = addfile(dir+nvar(iv)+"."+ssn+".nc","r")
  a = addfile("../output/matlab/yearly1/vars200.1980cl.nc","r")
  ;************************************************
  ;if(any((/"u","v"/).eq.nvar(iv))) then
  ;  lev = (/-60, 60, 5/)
  ;end if
  ;if(any((/"um","vm"/).eq.nvar(iv))) then
;    lev = (/-120, 120, 10/)
;  end if
;  if(any((/"umy"/).eq.nvar(iv))) then
;    lev = (/-50, 50, 5/)
;  end if
;  if(any((/"umx","vmx","vmy"/).eq.nvar(iv))) then
;    lev = (/-20, 20, 2/)
;  end if
;  if(any((/"q"/).eq.nvar(iv))) then
;    lev = (/-5, 5, 0.25/)
;  end if
;  if(any((/"qx"/).eq.nvar(iv))) then
;    lev = (/-20, 20, 5/)
;  end if
;  if(any((/"qy","BetaM"/).eq.nvar(iv))) then
;    lev = (/-20, 20, 2/)
;  end if
;  if(any((/"qxx","qyy","qxy"/).eq.nvar(iv))) then
;    lev = (/-100, 100, 50/)
;  end if
lev = (/-30, 30, 5/)
  ;************************************************
  ;  var = a->$nvar(iv)$
    var = a->RWS
  ;************************************************
  ; create default plot
  ;************************************************
    wks = gsn_open_wks("png","../output/matlab/yearly1/RWS200cl")                  ; send graphics to PNG file

    res = True
    res@mpMinLatF            = -90      ; range to zoom in on
    res@mpMaxLatF            =  90.
    res@mpMinLonF            =  0.
    res@mpMaxLonF            = 360.

    res@cnFillOn             = True                  ; turn on color fill
    res@cnFillPalette        = "BlWhRe"              ; choose colormap



    res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
    res@cnMinLevelValF  = lev(0)			; set the minimum contour level
    res@cnMaxLevelValF  = lev(1)			; set the maximum contour level
    res@cnLevelSpacingF = lev(2)			; set the interval between contours

    res@lbLabelStride         = 2

    ;print(nvar(iv)+"")
    printVarSummary(var)

    print("plotting")
    res@gsnRightString = ""
    plot = gsn_csm_contour_map(wks,var(1:240,:),res)         ; create a default plot
    print("end plotting")

;end do
end
