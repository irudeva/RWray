begin
;************************************************
; read in netCDF file
;************************************************

 nvar = (/"Ks"/)
 ssn = (/"Sep"/)
; do iv=0,dimsizes(nvar)-1
 do iv=0,0
  a = addfile("../output/Ks/"+nvar(iv)+"."+ssn+".erain.nc","r")
  ;************************************************
  if(any((/"Ks"/).eq.nvar(iv))) then
    lev = (/-20, 20, 2/)
  end if
  ;************************************************
    var = a->$nvar(iv)$

    if (any(isnan_ieee(var))) then
        value = 1.e20
        replace_ieeenan (var, value, 0)
        var@_FillValue = value
    end if

    ;print(var(:10,100))
  ;************************************************
  ; create default plot
  ;************************************************
    wks = gsn_open_wks("png","../output/Ks/"+nvar(iv)+"."+ssn+".erain")                  ; send graphics to PNG file

    res = True
    res@mpMinLatF            = -90      ; range to zoom in on
    res@mpMaxLatF            =  90.
    res@mpMinLonF            =  0.
    res@mpMaxLonF            = 360.
    res@mpCenterLonF         = 180.

    res@cnFillOn             = True                  ; turn on color fill
    res@cnFillPalette        = "BlWhRe"              ; choose colormap


    res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
    res@cnMinLevelValF  = lev(0)			; set the minimum contour level
    res@cnMaxLevelValF  = lev(1)			; set the maximum contour level
    res@cnLevelSpacingF = lev(2)			; set the interval between contours

    res@lbLabelStride         = 2

    print(nvar(iv)+"")
    printVarSummary(var)

    print("plotting")
    plot = gsn_csm_contour_map(wks,var(5:236,:),res)         ; create a default plot
    print("end plotting")

end do
end
