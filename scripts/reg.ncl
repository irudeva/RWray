;**********************************************************************
; polyg_13.ncl
;
; Concepts illustrated:
;   - Drawing filled boxes on a map plot using different methods
;   - Attaching filled boxes to a map plot
;   - Zooming in on a particular area on a Lambert Conformal map
;   - Turning on edges for polygons
;   - Attaching filled polygons to a map
;   - Changing the color of a filled polygon
;   - Using gc_latlon to create a straight box on a map
;   - Turning off the map lat/lon grid lines
;
;**********************************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;
;**********************************************************************
; This function attaches a filled box to a map.
;
; It gets the edges of the box using gc_latlon to calculate the
; great circle between two lat/lon points.
;
; This ensures that your box edges will be straight when you attach
; them to the map.
;
;**********************************************************************

function add_filled_box(wks,map,lftlat,lftlon,rgtlat,rgtlon, color)
local gnres, npts, i, lat_begend, lon_begend, dist
begin
  gnres             = True
  gnres@gsFillColor = color
  gnres@gsEdgesOn   = True

  npts   =  5     ; Number of points along each box edge.
                  ; You could make this different for each
                  ; edge, if you want.
;
; Define the coordinates for the start, end of the four sides
; of each box. Put them in a big array so it's easier to
; loop across the points later.
;
                         ; bottom            right
                         ; top               left
  lat_begend = (/ (/lftlat,lftlat/), (/lftlat,rgtlat/),  \
                  (/rgtlat,rgtlat/), (/rgtlat,lftlat/)/)
  lon_begend = (/ (/lftlon,rgtlon/), (/rgtlon,rgtlon/),  \
                  (/rgtlon,lftlon/), (/lftlon,lftlon/)/)

;---Define array to hold box.
  latbox = new(4*npts,float)
  lonbox = new(4*npts,float)

;---Loop across the four edges and calculate the points along each edge
  do i=0,3
    ibeg = i*npts
    iend = ibeg+npts-1
    dist = gc_latlon(lat_begend(i,0),lon_begend(i,0), \
                     lat_begend(i,1),lon_begend(i,1),npts,2)
    latbox(ibeg:iend) = dist@gclat
    lonbox(ibeg:iend) = dist@gclon
  end do

;---Attach filled box to map and return
  dum = gsn_add_polygon(wks, map, lonbox, latbox, gnres)

  return(dum)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
lonreg1 = (/ 5, 30, 70, 100,  80, 100, 230, 260, 30 ,160, 180 /)
lonreg2 = (/25, 50, 90, 120, 100, 120, 255, 280, 70, 180, 200 /)
latreg1 = (/40, 45, 50,  55,  25,  22,  45,  35, 70,  70,  70 /)
latreg2 = (/55, 60, 65,  70,  45,  42,  60,  50, 80,  80,  80 /)

nreg = dimsizes(lonreg1)

color=(/"yellow", "darkolivegreen1","coral2"/)

SIEcorr_JFM       = (/-.19, -0.16, -.04, -.02, -.39, -0.37, -0.15, -0.16, -0.70, -.66/)
SIEcorrsign_JFM   = (/   0,     0,    0,    0,    1,     1,     0,     0,     1,    1/)
SIEcorrdt_JFM     = (/ .22,  0.14,  .04,  .02,  .06,  0.27, -0.16, -0.07, -0.53, -.48/)
SIEcorrdtsign_JFM = (/   0,     0,    0,    0,    0,     0,     0,     0,     1,    1/)

SIEcorr_JAS       = (/-.51, -0.66, -.33, -.24, -.74, -0.64, -0.42, -0.41, -0.45, -.79/)
SIEcorrsign_JAS   = (/   1,     1,    1,    0,    1,     1,     1,     1,     1,    1/)
SIEcorrdt_JAS     = (/ .05, -0.17, -.43, -.15, -.08,  0.06, -0.05, -0.39, -0.33, -.50/)
SIEcorrdtsign_JAS = (/   0,     0,    1,    0,    0,     0,     0,     1,     1,    1/)


; dirout = "../output/SIE/"
dirout = "../output/regmap/"

dum = new(nreg,graphic)
dumcorr = new((/2,nreg/),graphic)



  wks_type = "png"
  wks_type@wkWidth = 2500
  wks_type@wkHeight = 2500

  wks = gsn_open_wks(wks_type,dirout+"regmap_ESib")  ; send graphics to PNG file

  res = True

  res@gsnDraw               = False
  res@gsnFrame              = False
  res@gsnMaximize           = True

; Add map resources
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  ;res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks

  res@mpProjection      = "Robinson"       ; choose projection
  res@mpGeophysicalLineThicknessF = 3.
  ;res@mpMinLatF   = 0.
  ;res@mpMaxLatF   = 90.

  res@mpLimitMode        = "Corners"
    res@mpLeftCornerLatF   = 0.
  ;  res@mpLeftCornerLonF   = -170.
    res@mpRightCornerLatF  = 84.
  ;  res@mpRightCornerLonF  = 170.

  res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpGridLatSpacingF = 30.               ; spacing for lat lines
  res@mpGridLonSpacingF = 30.               ; spacing for lon lines
  res@mpFillOn          = False


  ;res@mpLambertParallel1F = 30.
  ;res@mpLambertParallel2F = 60.
  ;res@mpLambertMeridianF  = 30.

;---Draw a box using gsn_add_polygon. See how the edges are curved.
  gnres             = True
  gnres@gsFillColor = "yellow"
  gnres@gsEdgesOn   = True

;---Create new map plot with a new title.

  res@tiMainString = "Regions and correlation of detrended t2m with raw SIE"
  map = gsn_csm_map(wks, res)

  tres                      = True                ; text mods desired
  tres@txFontHeightF        = 0.012               ; make smaller


;
; Calculate box edges using gc_latlon. Now we can attach box to
; the map, and the edges will be straight.
;
do ir = 0,nreg-4
  lftlat = latreg1(ir)
  rgtlat = latreg2(ir)
  lftlon = lonreg1(ir)
  rgtlon = lonreg2(ir)

  if(ir.le.nreg-4)then
    ic =0
  else
    ic =1
  end if
  if(ir.eq.3)then
     ic = 2
  end if

  dum(ir) = add_filled_box(wks,map,lftlat,lftlon,rgtlat,rgtlon, color(ic))

  if(ir.ne.nreg-1)then
  tres@txFontColor = "steelblue2"
  if(SIEcorrdtsign_JFM(ir).eq.1)then
    tres@txFontColor = "red"
  end if
;  dumcorr(0,ir) = gsn_add_text(wks,map,sprintf("%3.2f", SIEcorrdt_JFM(ir)),\
;                          lonreg1(ir)+10.,lftlat+12.,tres)
  tres@txFontColor = "steelblue2"
  if(SIEcorrdtsign_JAS(ir).eq.1)then
    tres@txFontColor = "red"
  end if
;  dumcorr(1,ir) = gsn_add_text(wks,map,sprintf("%3.2f", SIEcorrdt_JAS(ir)),\
;                          lonreg1(ir)+10.,lftlat+5.,tres)
  end if
end do

  draw(map)
  frame(wks)

end
