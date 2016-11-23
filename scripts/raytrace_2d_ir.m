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

rad     = 6.371e6  ; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5 ; % rotation rate of Earth (rad/s)
dtr     = pi/180   ;
rtd     = 180/pi   ;

%% specify the number of points at the north/south pole to remove from analysis
j_pole=5;



day = 24*60*60; % in s
dt  = 15*60;  % in s
integration_time=10*day;
Nsteps = round(integration_time/dt);

% Periods ('Inf' for quasi-stationary waves) 
%Periods=[ Inf 50 20 ]*day;
Periods=[ Inf ]*day;
freq=2*pi./Periods;
Nfr=length(freq);

% Initial k wave number:
k_wavenumbers=[3];
Nk = length(k_wavenumbers);

% Starting point of ray (location)

lon0 = [120 30 90]  ; %deg.E
lat0 = [50 30 5 ] ; %deg.N

% smoothing before ray tracing MIGHT be a good idea...:
do_smooth_background_fields=1;

do_complex_tracing=0

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
lev       = find(level==300);
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
    nlat    = size(lat,1)
    jmin    = 1 + j_pole
    jmax    = nlat - j_pole 
     display(['size of lat = ',sprintf(' %d',nlat)]);
else
    disp('Check lat');
    exit
end;
[name,type,dimids,natts] = netcdf.inqVar(ncid,0);
if(name == 'longitude'); 
    lon     = netcdf.getVar (ncid,0);
    lonin   = lon;
    nlon    = size(lon,1);
     display(['size of lon = ',sprintf(' %d',nlon)]);
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
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    psi = single(psi)*scale+offset;
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
for yr = 1980:2010
    for imon = 1:3
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

xm=transpose(rad*lon*dtr);
ym=rad*log((1+sin(dtr*lat))./cos(dtr*lat));
ym(lat==90)=inf; % was using 1e10 with success here
ym(lat==-90)=-inf;

%  Smoothing the fields if do_smooth_background_fields=1
if do_smooth_background_fields
   u0=zfltr(u0,1,10,1);
   v0=zfltr(v0,1,10,1);
   psi0=zfltr(psi0,1,10,1);
  %u00=smooth(u0,'loess');
  %u01=reshape(u00,nlat,nlon);
  %v00=smooth(v0,'loess');
  %v01=reshape(v00,nlat,nlon);
  %psi00=smooth(psi0,'loess');
  %psi01=reshape(psi00,nlat,nlon);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for UbarM, VbarM, and even psiM

[X,Y]=meshgrid(lon,lat);
[Lon,Lat]=meshgrid(lon*dtr,lat*dtr);

UbarM=u0./cos(Lat);
VbarM=v0./cos(Lat);
%while testing only PsiM=psi0;
PsiM=psi0./cos(Lat);

UbarM(Y==90)=inf; % was using 1e10 with success here
UbarM(Y==-90)=-inf;
VbarM(Y==90)=inf; % was using 1e10 with success here
VbarM(Y==-90)=-inf;
PsiM(Y==90)=inf; % was using 1e10 with success here
PsiM(Y==-90)=-inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for BetaM; NOTE that cos2(Lat)=(1+cos(2*Lat))/2
% 1)  BetaM--this can be solved for a number of ways, as Karoly did 
%%  (and also Hoskins and Karoly, 1981), where BetaM=2*Omega*cos2(Lat)/r
%%         - d/dy[(1/cos2(Lat))*d/dy(cos2(Lat)*UbarM)]
%%     or as Hoskins and Ambrizzi specify (2 ways: w/ and w/o Mercator) 

b=(1+cos(2*Lat))/2;
trm1=2*e_omega*b/rad;
c=b.*UbarM;
d=1./b;
cdy=c;cdy2=c;
for j=2:nlat-1
  cdy(j,:)=(mean(c(j-1:j,:))-mean(c(j:j+1,:)))/(mean(ym(j-1:j))- ...
                                                mean(ym(j:j+1)));
end
cdy=cdy.*d;
for j=3:nlat-2
  cdy2(j,:)=(mean(cdy(j-1:j,:))-mean(cdy(j:j+1,:)))/(mean(ym(j- ...
                                                    1:j))-mean(ym(j:j+1)));
end

BetaM=trm1-cdy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GRADIENTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for qbar
%%%%  The larger structure is to make sure have entire zonal field
%%%%  qbar = (d2/dx+d2/dy)psi/cos2(lat) + 2*Omega*sin(lat) - absolute vort on the sphere

tempPSI=[PsiM(:,nlon-3:nlon) PsiM PsiM(:,1:4)];
tempxm=[xm(nlon-3:nlon) xm xm(1:4)];

px1=NaN(size(tempPSI), 'like',tempPSI);px2=px1;
py1=px1;py2=px1;
for j=2:nlat-1
  py1(j,:)=(mean(tempPSI(j-1:j,:))-mean(tempPSI(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end
for j=3:nlat-2
  py2(j,:)=(mean(py1(j-1:j,:))-mean(py1(j:j+1,:)))/(mean(ym(j-1: ...
                                                    j))-mean(ym(j:j+1)));
end
for i=2:nlon+7
  px1(:,i)=(mean(tempPSI(:,i-1:i),2)-mean(tempPSI(:,i:i+1),2))/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
end
for i=3:nlon+6
  px2(:,i)=(mean(px1(:,i-1:i),2)-mean(px1(:,i:i+1),2))/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
end

dd=[d(:,nlon-3:nlon) d d(:,1:4)];  % d=1/cos2(lat)
f=(2*e_omega)*sin(Lat);
tempf=[f(:,nlon-3:nlon) f f(:,1:4)];

tempqbar=tempf+dd.*(px2+py2);
qbar=tempqbar(:,5:nlon+4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dqbar/dx, and d2qbar/dx2

clear px1 px2
px1=NaN*ones(nlat,nlon+8); px2=px1;
for i=4:nlon+5
  px1(:,i)=(mean(tempqbar(:,i-1:i),2)-mean(tempqbar(:,i:i+ ...
                                                    1),2))/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
  %display(['i=',sprintf(' %d',i)])
  %display(['xm=',sprintf(' %d',mean(tempxm(i-1:i))-mean(tempxm(i:i+1)))])
  %display(['tempqbar=',sprintf(' %d',tempqbar(10,i-1:i+1))])
  %display(['tempqbar1=',sprintf(' %d',mean(tempqbar(10,i-1:i),2))])
  %display(['tempqbar2=',sprintf(' %d',mean(tempqbar(10,i:i+1),2))])
  %display(['px1=',sprintf(' %d',px1(10,i))])
  
end
for i=5:nlon+4
  px2(:,i)=(mean(px1(:,i-1:i)')'-mean(px1(:,i:i+1)')')/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
end

dqbardx=px1(:,5:nlon+4);
d2qbardx2=px2(:,5:nlon+4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dqbar/dy (should be same as BetaM), and d2qbar/dy2

py1=NaN*ones(nlat,nlon+8); py2=py1;
for j=(jmin-1):(jmax+1)
  py1(j,:)=(mean(tempqbar(j-1:j,:))-mean(tempqbar(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end
for j=jmin:jmax
  py2(j,:)=(mean(py1(j-1:j,:))-mean(py1(j:j+1,:)))/(mean(ym(j-1: ...
                                                    j))-mean(ym(j:j+1)));
end

dqbardy=py1(:,5:nlon+4);
d2qbardy2a=py2(:,5:nlon+4);

%%%%  Alternately for d2qbar/dy2

py2=NaN*ones(nlat,nlon);
for j=4:nlat-3
  py2(j,:)=(mean(BetaM(j-1:j,:))-mean(BetaM(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end

d2qbardy2=py2;

%% A debugging text for the calculation of d^2 qbar/dy^2:
% sum(sum(d2qbardy2(6:68,6:138)-d2qbardy2a(6:68,6:138)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for d2qbar/dxy (once);  COULD ALSO DO dx of either BetaM or
%%%%  dqbardy

py1=NaN*ones(nlat,nlon);
for j=2:nlat-1
  py1(j,:)=(mean(dqbardx(j-1:j,:))-mean(dqbardx(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end

d2qbardxy=py1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dUbarM/dx and dUbarM/dy

tempU=[UbarM(:,nlon) UbarM UbarM(:,1)];

px1=NaN*ones(nlat,nlon+2);py1=px1;
for j=2:nlat-1
  py1(j,:)=(mean(tempU(j-1:j,:))-mean(tempU(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end
for i=2:nlon+1
  px1(:,i)=(mean(tempU(:,i-1:i)')'-mean(tempU(:,i:i+1)')')/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
end

dUbarMdx=px1(:,2:nlon+1);
dUbarMdy=py1(:,2:nlon+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dVbarM/dx and dVbarM/dy

tempV=[VbarM(:,nlon) VbarM VbarM(:,1)];

px1=NaN*ones(nlat,nlon+2);py1=px1;
for j=2:nlat-1
  py1(j,:)=(mean(tempV(j-1:j,:))-mean(tempV(j:j+1,:)))/ ...
           (mean(ym(j-1:j))-mean(ym(j:j+1)));
end
for i=2:nlon+1
  px1(:,i)=(mean(tempV(:,i-1:i)')'-mean(tempV(:,i:i+1)')')/ ...
           (mean(tempxm(i-1:i))-mean(tempxm(i:i+1)));
end

dVbarMdx=px1(:,2:nlon+1);
dVbarMdy=py1(:,2:nlon+1);

%%%%%%%%%%%%%%%%%%%%%
%% Start ray tracing:
%%%%%%%%%%%%%%%%%%%%%

%% Solving for the ray path for different forcing sites (initial
%% locations of rays):

for iloc=2:2 %size(lon0,2)

    [tmp,i0] = min(abs(lon-lon0(iloc)));
    [tmp,j0] = min(abs(lat-lat0(iloc)));
    
    
    [sxm,sym]=meshgrid(xm,ym(jmin:jmax));
    [XX,YY]=meshgrid(lon,lat(jmin:jmax));
    sBetaM=BetaM(jmin:jmax,:);
    sUbarM=UbarM(jmin:jmax,:);
    sVbarM=VbarM(jmin:jmax,:);
    sqbar=qbar(jmin:jmax,:);
    sdqbardy=dqbardy(jmin:jmax,:);
    sdqbardx=dqbardx(jmin:jmax,:);
    sd2qbardx2=d2qbardx2(jmin:jmax,:);
    sd2qbardy2=d2qbardy2(jmin:jmax,:);
    sd2qbardxy=d2qbardxy(jmin:jmax,:);
    sdUbarMdx=dUbarMdx(jmin:jmax,:);
    sdUbarMdy=dUbarMdy(jmin:jmax,:);
    sdVbarMdx=dVbarMdx(jmin:jmax,:);
    sdVbarMdy=dVbarMdy(jmin:jmax,:);
  
    
    %%%%  Estimating the initial Ks from the forcing site
  
    for iomega=1:Nfr
     omega=freq(iomega);
     period=round((2*pi/omega)/day);
     
     
    
     for kkr=1:Nk
      kk=k_wavenumbers(kkr);
      fprintf('Ray tracing... period=%d, k=%d, ilocation=%d\n' ...
              ,period,kk,iloc)
      spotk=kk/rad/cos(Lat(j0,i0));
      %subk=[ilocation kk RR];
      yi=ym(j0);xi=xm(i0);
      %%  Ks=sqrt(abs(sBetaM(fry,frx)/sUbarM(fry,frx)));
      %%  spotl=real((Ks^2-spotk^2)^0.5);
      Bint=BetaM(j0,i0);Uint=UbarM(j0,i0);Vint=VbarM(j0,i0);
      qxint=dqbardx(j0,i0);qyint=dqbardy(j0,i0);
      %% Calculate the initial l wave number from the initial omega
      %% and k by solving the polynomial equation based on the
      %% dispersion relation (equation 8 in Karoly 1983):
      %% change the following to have a non zero frequency:
        
      %% Eli: added omega terms here for the non zero frequency case:
      cz(1)=Vint;
      cz(2)=Uint*spotk-omega;
      cz(3)=Vint*spotk^2+qxint;
      cz(4)=Uint*spotk^3-qyint*spotk-omega*spotk^2;
      tl=roots(cz);
      %tl*rad
      %Vint
      %Uint
      %omega
      %spotk*rad*cos(Lat(j0,1))
      %qxint
      %qyint
        
      for RR=1:3  
        spotl=tl(RR);
        Ks=(spotl^2+spotk^2)^0.5;
 
        fprintf('Ray tracing... period=%d, k=%d, root=%d, l=%d, ilocation=%d\n' ...
                ,period,kk,RR,spotl*rad,iloc);

        %if do_only_northern_hisphere_rays
          %%  If only interested in rays that go northward, break the
          %%  run if the starting l wavenumber has negative real part:
          %tstl=real(spotl)*r*cos(Lat(fry+4,1));
          %if real(spotl)<=0 
            %fprintf(1,['*** found real(l)<0, ray going to southern ' ...
            %           'hemisphere, not tracing it\n']);
            %fprintf('[real(l) ilocation kk RR]=[%g %g %g %g]' ...
            %        ,real(spotl),ilocation,kk,RR);
            %break 
          %end
        %end
        
        if do_complex_tracing==0
          %%  If not interested in complex ray tracing, break the run
          %%  if the starting l wavenumber has an imaginary part.  N.B.
          %%  The complex conjugates give redundant solutions, so only
          %%  need to trace one
          tstl2=imag(spotl)*rad*cos(Lat(j,1));
          if tstl2~=0
            fprintf('*** found complex initial l, not tracing. \n')
            fprintf('    [tstl2 ilocation kk omega RR]=[%g,%g,%g,%g,%g]\n' ...
                    ,tstl2,iloc,kk,omega,RR)
            continue
          end
        end
        
        %%%%% output file %%%%%%%%%%%%%
        %fn_out = sprintf(...
        %     '../output/matlab/ray_lat%dN__lon%dN_period%d_k%d_root%d',...
        %     lat0(iloc),lon0(iloc),period,kk,RR);
  
        %%%%%%%%%%%
        %%  Starting the loop with the above initial k,l, and Ks
        
        for t=1:Nsteps
          if rem(t,day/(2*dt))==0
            fprintf(1,'t = %g\n',t);
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%  Interpolating the fields to the current spot      

          Uint=interp2(sxm,sym,sUbarM,xi,yi,'spline');
          Vint=interp2(sxm,sym,sVbarM,xi,yi,'spline');
          Uxint=interp2(sxm,sym,sdUbarMdx,xi,yi,'spline');
          Uyint=interp2(sxm,sym,sdUbarMdy,xi,yi,'spline');
          Vxint=interp2(sxm,sym,sdVbarMdx,xi,yi,'spline');
          Vyint=interp2(sxm,sym,sdVbarMdy,xi,yi,'spline');
          qint=interp2(sxm,sym,sqbar,xi,yi,'spline');
          qxint=interp2(sxm,sym,sdqbardx,xi,yi,'spline');
          qyint=interp2(sxm,sym,sdqbardy,xi,yi,'spline');
          qxxint=interp2(sxm,sym,sd2qbardx2,xi,yi,'spline');
          qyyint=interp2(sxm,sym,sd2qbardy2,xi,yi,'spline');
          qxyint=interp2(sxm,sym,sd2qbardxy,xi,yi,'spline');
          Bint=interp2(sxm,sym,sBetaM,xi,yi,'spline');
          Cosint=interp2(sxm,sym,cos(Lat(jmin:jmax,:)),xi,yi,'spline');
           
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%  Solving for the changes
          
          dkdt=-spotk*Uxint-spotl*Vxint+(qxyint*spotk-qxxint*spotl)/Ks^2;
          dldt=-spotk*Uyint-spotl*Vyint+(qyyint*spotk-qxyint*spotl)/Ks^2;
          dxdt=Uint+((spotk^2-spotl^2)*qyint-2*spotk*spotl*qxint)/Ks^4;
          dydt=Vint+(2*spotk*spotl*qyint-(spotk^2-spotl^2)*qxint)/Ks^4;
          
          %%%%%%%%%%%%%%
          %%  Updating the changes
          
          xi=xi+real(dxdt)*dt;
          if xi>=max(max(sxm))
            xi=xi-max(max(sxm));
          end
          yi=yi+real(dydt)*dt;
          spotl=spotl+dldt*dt;
          spotk=spotk+dkdt*dt;
          Ks=(spotk^2+spotl^2)^0.5;
          
          if rem(t,day/(2*dt))==0
            disp(['xi=',sprintf(' %d',xi),' yi=',sprintf(' %d',yi),...
              ' spotk=',sprintf(' %d',spotk*rad*Cosint),' spotl=',sprintf(' %d',spotl*rad)])
          end
          
          %%%%%%%%%%%%%%%
          %%  Finding the location
          
           Yint=interp2(sxm,sym,YY,xi,yi,'spline');
           Xint=interp2(sxm,sym,XX,xi,yi,'spline');
 
          %% make sure ray does not leave the domain where
          %% background fields are given:
          if Yint<lat(jmax) || Yint>lat(jmin)
            fprintf(1,'*** Yint>Ymax, breaking.  (Xint,Yint)=(%g,%g)\n'...
                    ,Xint,Yint);
            break
          end
          
         %%%%%%%%%%%%%%
          %%  Storing
          adj=rad*Cosint;
          timestep(t,:) =[t t*dt/3600 t*dt/day];
          locm(t,:)=[xi yi];
          rwnums(t,:)=[real(spotk)*adj real(spotl)*rad];
          iwnums(t,:)=[imag(spotk)*adj imag(spotl)*rad];
          rw(t,:)=[real(dxdt) real(dydt)];
          iw(t,:)=[imag(dxdt) imag(dydt)];
          locgeo(t,:)=[Xint Yint];
          %wchg(t,:)=[real(dldt)*adj real(dkdt)*adj imag(dldt)*adj imag(dkdt)*adj];
          %omega
          %rsom(t,:)=[real(Uint*spotk+Vint*spotl+(qxint*spotl-qyint*spotk)/Ks^2)];
          %isom(t,:)=[imag(Uint*spotk+Vint*spotl+(qxint*spotl-qyint*spotk)/Ks^2)];
        end
        for t = 1:Nsteps
          if rem(t,day/(2*dt))==0
              %disp([time, locm, locgeo, rwnums, iwnums])
            alL(2*t*dt/day,:)=[timestep(t,:) locm(t,:) locgeo(t,:) rwnums(t,:) iwnums(t,:)];
          end
         end
         fn_out = sprintf('../output/matlab/ray_lat%dN__lon%dN_period%d_k%d_root%d',...
                         lat0(iloc),lon0(iloc),period,kk,RR);
         dlmwrite(fn_out, alL,'precision', '%.6f');
            
      end
     end
    end
end



%------------------------------------------------
%%%  !!! add time from the start of inegration (s) to the output !!!
%fout = sprintf('../output/ray_loc%dE_%dN_period%d_k%d_root%d'...
%                         ,lon0(),lat0(),period,kk,RR);
