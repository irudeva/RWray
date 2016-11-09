%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate 2d Barotropic Rossby ray paths over specified background 
%%  fields.
%
%   Based on Karoly, 1983
%   Varables:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear all

addpath matlab matlab/spherepack

k=1

rad     = 6.371e6  ; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5 ; % rotation rate of Earth (rad/s)
dtr     = pi/180   ;
rtd     = 180/pi   ;

day = 24*60*60; % in s
dt  = 60*60;  %
integration_time=10*day;
Nsteps = round(integration_time/dt);

% Periods ('Inf' for quasi-stationary waves) 
Periods=[ Inf 50 20 ]*day;
freq=2*pi./Periods;
Nfrequencies=length(freq);

% Initial k wave number:
k_wavenumbers=[1 2 3 4 5 6];

% Starting point of ray (location)

lon0 = [ 0   0]  ; %deg.E
lat0 = [30 -30] ; %deg.N

% smoothing before ray tracing MIGHT be a good idea...:
do_smooth_background_fields=0;

disp('-----------------------------------');

disp(['k wavenumbers: ' num2str(k_wavenumbers(:)') '']) ;
disp(['Periods: ' num2str(Periods(:)'/day) ' days']) ;
disp(['Integration time: ', num2str(integration_time'/day),' days']) ;
disp(['dt: ' num2str(dt'/60) ' min']) ;
disp(['Stating point: ' num2str(lon0) ' E   ' num2str(lat0) ' N']) ;

disp('-----------------------------------');
disp('');
% Read data

fuwnd     = '../data/wnd300.mnth.erain.nc';
disp(['U wind read from ' fuwnd])
ncid      = netcdf.open ( fuwnd,'NC_NOWRITE' );
level     = netcdf.getVar (ncid,0);
lev       = find(level==200);
varid     = 3;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
if(name=='u'); 
    uwnd     = netcdf.getVar (ncid,varid,'short');
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    uwnd = single(uwnd)*scale+offset;   
    %u = squeeze(uwnd(:,:,lev,:));
    u = uwnd;
     display(['size of uwnd (',name,') = ',sprintf(' %d',size(u))]);
end;
[name,type,dimids,natts] = netcdf.inqVar(ncid,1);
if(name == 'latitude'); 
    lat     = netcdf.getVar (ncid,1);
    latin   = lat;
     display(['size of lat = ',sprintf(' %d',size(lat))]);
else
    disp('Check lat');
    exit
end;
[name,type,dimids,natts] = netcdf.inqVar(ncid,0);
if(name == 'longitude'); 
    lon     = netcdf.getVar (ncid,0);
    lonin   = lon;
     display(['size of lon = ',sprintf(' %d',size(lon))]);
else
    disp('Check lon');
    exit
end;
[name,type,dimids,natts] = netcdf.inqVar(ncid,2);
if(name == 'time'); 
    time     = netcdf.getVar (ncid,2);
     display(['size of time = ',sprintf(' %d',size(time))]);
else
    disp('Check time');
    exit
end;



fvwnd     = '../data/wnd300.mnth.erain.nc';
disp(['V wind read from ' fvwnd])
ncid      = netcdf.open ( fvwnd,'NC_NOWRITE' );
varid     = 4;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
if(name=='v'); 

    vwnd     = netcdf.getVar (ncid,varid,'short');
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    vwnd = single(vwnd)*scale+offset;
    %v = squeeze(vwnd(:,:,lev,:));
    v = vwnd; 
     display(['size of vwnd (',name,') = ',sprintf(' %d',size(v))]);
end;

fpsi      = '../data/sf300.mnth.erain.nc' ;
disp(['Streamfunction read from ' fpsi])
ncid      = netcdf.open ( fpsi,'NC_NOWRITE' );
varid     = 3;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
if(name=='sf'); 
    %disp('reading streamfunction');
    psi     = netcdf.getVar (ncid,varid);  
     display(['size of streamfunction (',name,') = ',sprintf(' %d',size(psi))]);
end;
netcdf.close(ncid);

disp('done getting data.');

% Select time for the background field (climatology)

TIMEbase = datenum(1900, 1, 1);
date = datestr(double(time/24) + TIMEbase); % where time is hours from Timebase

formatdata = '01-%s-%d';
mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!
in=1;
for yr = 1980:1980
    for imon = 1:1
      nyr = yr;
      if(mon(imon,:)=='Dec'); nyr=yr-1; end
      mydate(in,:) = sprintf(formatdata,mon(imon,:),nyr);
      %display(['My    date: ',mydate(in,:)])
      in = in+1;
    end
end
ndate = size(mydate,1);
 
for in = 1:ndate
 for t = 1:size(time,1)
  date = datestr(double(time(t)/24) + TIMEbase);
   if date == mydate(in,:)
    display(['Found date: ',date])
    nt(in)  = t;
  end
 end
end

% Climatology
u0=(squeeze(mean(u(:,:,nt),3)))';
[m,n]=size(u0);
v0=(squeeze(mean(v(:,:,nt),3)))';
psi0=(squeeze(mean(psi(:,:,nt),3)))';

% Check that data are from N to S
if lat(1)<lat(end)
    if size(u0,1)~=length(lat) || size(u0,2)~=length(lon)
      disp('Flipud: u0 array mismatched to dimensions')
      exit
    end
    lat=flipud(lat);
    u0=flipud(u0);
    v0=flipud(v0);
    psi0=flipud(psi0);
    disp('The data were flipped to be N->S')
end


% In Mercator projection

xm=rad*lon*dtr;
ym=rad*log((1+sin(dtr*lat))./cos(dtr*lat));
ym(lat==90)=inf; % was using 1e10 with success here
ym(lat==-90)=-inf;

%  Smoothing the fields if do_smooth_background_fields=1
if do_smooth_background_fields
  u0=zfltr(u0,1,10,1);
  v0=zfltr(v0,1,10,1);
  psi0=zfltr(psi0,1,10,1);
end

%% Calculate all required gradients on Mercator projection

[ux,uy]=gradsphere(lon,lat,u0);



%------------------------------------------------
%%%  !!! add time from the start of inegration (s) to the output !!!
%fout = sprintf('../output/ray_loc%dE_%dN_period%d_k%d_root%d'...
%                         ,lon0(),lat0(),period,kk,RR);
