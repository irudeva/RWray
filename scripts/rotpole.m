function [nlat, nlon]  = rotpole(plat,plon,ilat,ilon)
% the location of the new pole is (plat,plon)
% the new coordinates are in a Cartesian c. system where:
% real North pole is (y=90,x=0) in the new system if plat > 0 and the East is (y=0,x=90)
% real South pole is (90,0) in the new system if plat < 0 and the East is (y=0,x=-90)

    rad  = 637103; %in m
    dtr  = pi/180.;
    rtd  = 180./pi;

    % copy initial arrays
    plat0 = plat;
    plon0 = plon;

    plat = (plat-90.)*dtr;
    plon = plon*dtr;

    ilon = ilon*dtr;
    ilat = ilat*dtr;

    nlat = zeros('like',ilat);
    nlon = zeros('like',ilon);
    

    for i=1:size(ilon,2)
    %for index, (lon,lat) in enumerate(zip(ilon, ilat)):
 
        % convert to cartesian coordinates
        % from bronstein p.217
        cos_plon = cos(plon);
        if plon == 90*dtr || plon == -90*dtr 
            cos_plon = 0;
        end
        cos_lon = cos(ilon(i));
        if ilon(i) == 90*dtr || ilon(i) == -90*dtr 
            cos_lon = 0;
        end

        x=rad*cos(ilat(i))*cos_lon;
        y=rad*cos(ilat(i))*sin(ilon(i));
        z=rad*sin(ilat(i));

        % turn XY
        xn=x*cos_plon+y*sin(plon);
        yn=x*(-sin(plon))+y*cos_plon;
        zn=z;

        % tilt XZ
        xnn=xn*cos(plat)+zn*sin(plat);
        ynn=yn;
        znn=xn*(-sin(plat))+zn*cos(plat);

        % convert back to polar coordinates
        rr=sqrt(xnn*xnn+ynn*ynn+znn*znn);
        nlat(i) = asin(znn/rr)*rtd;


        nlon(i) = atan2(xnn,ynn)*rtd;
        
     end
