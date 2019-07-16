%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate 2d Barotropic Rossby ray paths over specified background 
%%  fields.
%
%   Based on Karoly, 1983
%   Varables:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; 
% clear all

addpath matlab matlab/spherepack

rad     = 6.371e6  ; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5 ; % rotation rate of Earth (rad/s)
dtr     = pi/180   ;
rtd     = 180/pi   ;

%% specify the number of points at the north/south pole to remove from analysis
j_pole=5;

fyr = 1979;
lyr = 1985;

bgf = ['DJF'] % it depends on the input background fields
if bgf=='DJF'
 mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!
elseif bgf=='JJA'
 mon=[ 'Jun';'Jul';'Aug' ];  %%%!!!! for DEC -  year = year-1!!!
elseif bgf=='MAM'
 mon=[ 'Mar';'Apr';'May' ];  
elseif bgf=='SON'
 mon=[ 'Sep';'Oct';'Nov' ];  
end

% parameter setting: u wind
 % uwnd_dim = 0 to select 1 dim wind (zonally aver) 
 % uwnd_dim = 1 to select 2 dim wind (real wind)
uwnd_dim = 1;

day = 24*60*60; % in s
dt  = 15*60;  % in s
integration_time=10*day;
Nsteps = round(integration_time/dt);

% Periods ('Inf' for quasi-stationary waves) 
%Periods=[ Inf 50 20 ]*day;
Periods=[ Inf ]*day;

freq=2*pi./Periods;
Nfr=length(freq);

kout = 24 % number of points per day in the output file

% % Initial k wave number:
k_wavenumbers=[1 2 3 4 5 6];
%k_wavenumbers=[ 4 ];
Nk = length(k_wavenumbers);

% Starting point of ray (location)

%lat0 = [ 90:-90:10 ] ; %deg.N
%lat0 = [80:-10:-80] ; %deg.N
%lon0 = 150.*ones(size(lat0)) ; deg E

% lat0 = [80:-5:60] ; %deg.N
% lon0 = [60:90:360] ; %deg E
%lat0 = [-17 -51] ; %deg.N
%lon0 = [105 60] ; %deg E
 lat0 = [65] ; %deg.N
 lon0 = [10] ; %deg E

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

for iyr = fyr:lyr  
fuwnd     = sprintf('../data/erain.sf300.monmean.%d.nc',iyr);
disp(['U wind read from ' fuwnd])
ncid      = netcdf.open ( fuwnd,'NC_NOWRITE' );
%level     = netcdf.getVar (ncid,0);
%lev       = find(level==300);

varid     = 3;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
if(name=='u'); 
    uwnd     = netcdf.getVar (ncid,varid,'short');
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    uwnd = single(uwnd)*scale+offset;   
    %u = squeeze(uwnd(:,:,lev,:));
    u(:,:,:,iyr - fyr+1) = uwnd;
     display(['size of uwnd (',name,') = ',sprintf(' %d',size(u))]);
else
    disp('Check u');
    exit
end;
varid     = 4;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
if(name=='v'); 

    vwnd     = netcdf.getVar (ncid,varid,'short');
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    vwnd = single(vwnd)*scale+offset;
    %v = squeeze(vwnd(:,:,lev,:));
    v(:,:,:,iyr - fyr+1) = vwnd; 
     display(['size of vwnd (',name,') = ',sprintf(' %d',size(v))]);
else
    disp('Check v');
    exit
end;
varid     = 5;
[name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
 if(name=='sf'); 
    %disp('reading streamfunction');
    psi1     = netcdf.getVar (ncid,varid);  
    scale = netcdf.getAtt(ncid,varid,'scale_factor');
    offset= netcdf.getAtt(ncid,varid,'add_offset');
    psi(:,:,:,iyr - fyr+1) = single(psi1)*scale+offset;
     display(['size of streamfunction (',name,') = ',sprintf(' %d',size(psi))]);
else
    disp('Check sf');
    exit
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
    time(:,iyr - fyr+1)     = netcdf.getVar (ncid,2);
     display(['size of time = ',sprintf(' %d',size(time))]);
else
    disp('Check time');
    exit
end;
netcdf.close(ncid);

end 

disp('done getting data.');

if (uwnd_dim == 0);
    % Averaging of u - will write zonally averaged u to the grid points
    if size(u,1)==length(lon) || size(u,2)==length(lat)
      zUbar = repmat(mean(u,1),[size(lon),1]); 
    end
end
if (uwnd_dim == 1);
    % without averaging
    zUbar = u;
end


% -------------------------------------------------------------------

% Select time for the background field (climatology)

TIMEbase = datenum(1900, 1, 1);
date = datestr(double(time/24) + TIMEbase); % where time is hours from Timebase

formatdata = '01-%s-%d';
%mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!
in=1;
fyr1 = fyr
if bgf == 'DJF'
    fyr1 = fyr+1;
end
for yr = fyr1:lyr
    for imon = 1:3
      nyr = yr;
      if(mon(imon,:)=='Dec'); nyr=yr-1; end
      if (nyr>=fyr) && (nyr<=lyr)
       mydate(in,:) = sprintf(formatdata,mon(imon,:),nyr);
       %display(['My    date: ',mydate(in,:)])
       in = in+1;
      end
    end
end
ndate = size(mydate,1);
 
for in = 1:ndate
 for t1 = 1:size(time,1)
     for t2 = 1:size(time,2)
       date = datestr(double(time(t1,t2)/24) + TIMEbase);
     if date == mydate(in,:)
      display(['Found date: ',date])
%       nt1(in)  = t1;
%       nt2(in)  = t2;
      ussn(:,:,in) = zUbar(:,:,t1,t2); % for zonally symmetric case 
%       vssn(:,:,in) = v(:,:,t1,t2);
%       psissn(:,:,in) = psi(:,:,t1,t2);
     end
     end
 end
end

% Climatology
u0=(squeeze(mean(ussn(:,:,:),3)))';
[m,n]=size(u0);
% v0=(squeeze(mean(vssn(:,:,:),3)))';
% psi0=(squeeze(mean(psissn(:,:,:),3)))';

% Check that data are from N to S
if lat(1)<lat(end)
    if size(u0,1)~=length(lat) || size(u0,2)~=length(lon)
      disp('Flipud: u0 array mismatched to dimensions')
      exit
    end
    lat=flipud(lat);
    u0=flipud(u0);
%     v0=flipud(v0);
%     psi0=flipud(psi0);
    disp('The data were flipped to be N->S')
end


% In Mercator projection

xm=transpose(rad*lon*dtr);
ym=rad*log((1+sin(dtr*lat))./cos(dtr*lat));
ym(lat==90)=inf; % was using 1e10 with success here
ym(lat==-90)=-inf;

%  Smoothing the fields if do_smooth_background_fields=1
if do_smooth_background_fields
    %test
   %u0=zfltr(u0,1,10,1);
   %v0=zfltr(v0,1,10,1);
   %psi0=zfltr(psi0,1,10,1);
   %end test
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
% VbarM=v0./cos(Lat);
% PsiM=psi0;
%PsiM=psi0./cos(Lat);

UbarM(Y==90)=inf; % was using 1e10 with success here
UbarM(Y==-90)=-inf;
% VbarM(Y==90)=inf; % was using 1e10 with success here
% VbarM(Y==-90)=-inf;
% PsiM(Y==90)=inf; % was using 1e10 with success here
% PsiM(Y==-90)=-inf;

if (uwnd_dim == 1);
    UbarM_2d = UbarM;
end

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
%cdy=c;cdy2=c;
%cdy = NaN
for j=3:nlat-2
  cdy(j,:)=(mean(c(j-1:j,:))-mean(c(j:j+1,:)))/(mean(ym(j-1:j))- ...
                                                mean(ym(j:j+1)));
                                            continue
end
cdy(1:2,:) = nan ;
cdy(nlat-1:nlat,:) = nan ;
cdy=cdy.*d;
for j=4:nlat-3
  cdy2(j,:)=(mean(cdy(j-1:j,:))-mean(cdy(j:j+1,:)))/(mean(ym(j- ...
                                                    1:j))-mean(ym(j:j+1)));
  continue
end
cdy2(1:3,:) = nan ;
cdy2(nlat-2:nlat,:) = nan ;

BetaM=trm1-cdy2;

if uwnd_dim == 1
    BetaM_2d = BetaM;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%% Start ray tracing:
%%%%%%%%%%%%%%%%%%%%%

%%%%% Zonally symmetric case:
%%%%%   Ks2 = BetaM/(UbarM-Omega/k)
%%%%%   l = sqrt(Ks2 - k2)
%%%%%   ug = Omega/k + 2*BetaM*k2/Ks4
%%%%%   vg = 2*BetaM*k*l/Ks4


%% Solving for the ray path for different forcing sites (initial
%% locations of rays):
for ilat=1:size(lat0,2)
% for ilon=1:size(lon0,2)
ilon = ilat

    [tmp,i0] = min(abs(lon-lon0(ilon)));
    [tmp,j0] = min(abs(lat-lat0(ilat)));
    for iomega=1:Nfr
     omega=freq(iomega);
     period=round((2*pi/omega)/day);
     
     [sxm,sym]=meshgrid(xm,ym(jmin:jmax));
     [XX,YY]=meshgrid(lon,lat(jmin:jmax));
     sBetaM=BetaM_2d(jmin:jmax,:);
     sUbarM=UbarM_2d(jmin:jmax,:);
     scosLat=cos(Lat(jmin:jmax,:));
     
     %masking
     ind1 = find(sUbarM<5);
     ind2 = find(sBetaM<0);

    sBetaM1 = sBetaM;
    sBetaM1(ind1) = NaN;
    sBetaM1(ind2) = NaN;
    sUbarM1 = sUbarM;
    sUbarM1(ind1) = NaN;
    sUbarM1(ind2) = NaN;

     sBetaM=BetaM(jmin:jmax,:);
     sUbarM=UbarM(jmin:jmax,:);
    
     for kkr=1:Nk
      kk=k_wavenumbers(kkr);
      fprintf('Ray tracing... period=%d, k=%d, ilocation=(%dE, %dN)\n' ...
              ,period,kk,lon0(ilon),lat0(ilat))
      
      spotk0=rdivide(kk/rad,scosLat);
      %subk=[ilocation kk RR];
%       yi=ym(j0-jmin+1);xi=xm(i0);
      
      sKs2=abs(rdivide(sBetaM1,sUbarM1-rdivide(omega,spotk0)));
      ind3 = find(isnan(sKs2));
      sKs2(ind3) = 0;
      spotl1=real(sqrt(sKs2-spotk0.*spotk0));
      spotl2=-spotl1;
      
      sUbarM(ind3) = 0;
      sUbarM(ind3) = 0;
      
      %Bint=BetaM(j0,i0);Uint=UbarM(j0,i0);Vint=VbarM(j0,i0);
      
      for RR=1:2  
        yi=ym(j0);xi=xm(i0);
        if RR==1 
            spotl0 = spotl1;  
            coef = 1;
        elseif RR ==2
            spotl0 = spotl1;
            coef = -1;
        end
        
%         ug=rdivide(omega,spotk0)+rdivide(times(2*sBetaM,spotk0.^2),sKs2.^2);
%         vg=rdivide(times(2*sBetaM,spotk0.*spotl0),sKs2.^2);
        
        fprintf('Ray tracing... period=%d, k=%d, root=%d, l=%d, ilocation=%dE %dN\n' ...
                ,period,kk,RR,spotl0(j0-jmin+1,i0)*rad,lon0(ilon),lat0(ilat));
            
% a parameter to help with reflection
         refl = 0;    
            
        for t=1:Nsteps
          %  fprintf(1,'t = %g\n',t);
          if rem(t,day/(2*dt))==0
            fprintf(1,'t = %g\n',t);
          end

          spotl = interp2(sxm,sym,spotl0,xi,yi,'spline');
          %fprintf(1,'spotl = %g\n',spotl*rad);
          spotk = interp2(sxm,sym,spotk0,xi,yi,'spline');
          %fprintf(1,'spotk = %g without cos(Lat) adj!\n',spotk*rad);
          %%Ks2 = spotk^2+spotl^2;    
            
          Uint=interp2(sxm,sym,sUbarM,xi,yi,'spline');
          Bint=interp2(sxm,sym,sBetaM,xi,yi,'spline');
          Cosint=interp2(sxm,sym,scosLat,xi,yi,'spline');
          Ks2int=interp2(sxm,sym,sKs2,xi,yi,'spline');
          
%           spotl = real(sqrt(Ks2int-spotk*spotk));
          
          if abs(spotl*rad) <0.1
             fprintf(1,'*** spotl<0.1, breaking: spotl=%d (Xint,Yint)=(%g,%g)\n', spotl*rad,Xint,Yint);
             t = t-1;
             break
          end

          
          
%           dxdt = interp2(sxm,sym,ug,xi,yi,'spline');
%           dydt = interp2(sxm,sym,vg,xi,yi,'spline');

          if abs(spotl*rad) < 1. & refl == 0
%               if (locm(t-1,2) > locm(t-2,2) & locm(t,2)> locm(t-1,2)) || (locm(t-1,2) < locm(t-2,2) & locm(t,2) < locm(t-1,2))
%             fprintf(1,'*** spotl<0.1, breaking: spotl=%d (Xint,Yint)=(%g,%g)\n', spotl*rad,Xint,Yint);
%             break
%%                 spotl0 = -spotl0;
%               end
                coef = -coef;
                refl = 1;
                refl
          end
           if abs(spotl*rad) >= 1.1 & refl == 1    % .99*abs(spotl0(j0-jmin+1,i0)) & refl == 1
               refl = 0;
               refl
           end

          dxdt = omega/spotk +2*Bint*spotk^2/Ks2int^2;
          dydt = coef *2*Bint*spotk*spotl/Ks2int^2;
 
          %%%%%%%%%%%%%%
          %%  Updating the changes
          
          xi=xi+real(dxdt)*dt;
          if xi>=max(max(sxm))
            xi=xi-max(max(sxm));
          end
          yi=yi+real(dydt)*dt;
          
          if rem(t,day/(2*dt))==0
            disp(['xi=',sprintf(' %d',xi),' yi=',sprintf(' %d',yi),...
              ' spotk=',sprintf(' %d',spotk*rad*Cosint),' spotl=',sprintf(' %d',spotl*rad)])
          end
          
          %%%%%%%%%%%%%%%
          %%  Finding the geo location
          
           Yint=interp2(sxm,sym,YY,xi,yi,'spline');
           Xint=interp2(sxm,sym,XX,xi,yi,'spline');
          if rem(t,day/(2*dt))==0
           fprintf(1,'t=%d (Xint,Yint)=(%g,%g) (spotk,spotl)=(%g,%g) \n'...
                    ,t, Xint,Yint,spotk*rad*Cosint,spotl*rad); 
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
          %rsom(t,:)=[real(Uint*spotk-Bint*spotk/Ks2)];
          %isom(t,:)=[imag(Uint*spotk-Bint*spotk/Ks2)];
          %wave phase
          %rphase(t,:) = [(real(dxdt)*real(spotk) + real(dydt)*real(spotl) - rsom(t,:)*t)*dt];
          if t>1 
              rphase(t,:) = rphase(t-1,:) + [((real(dxdt)*real(spotk) + real(dydt)*real(spotl)))*dt];
          else
              rphase(t,:) = [((real(dxdt)*real(spotk) + real(dydt)*real(spotl)))*dt];
          end
          rphase1 = rphase(t,:);
          if rphase1 > 2*pi 
              while rphase1 > 2*pi
                  rphase1 = rphase1 - 2*pi;
              end
          else if rphase1 < 0
              while rphase1 < 0
                  rphase1 = rphase1 + 2*pi;
              end
          end
          end
          rphase2pi(t,:) = rphase1;
          
% not used
%           if t>1
%               if rphase1 >= 0
%                   if (rphase(t-1,:) <= pi/2 && rphase1 > pi/2) | (rphase(t-1,:) < 3*pi/2 & rphase1 > 3*pi/2)
%                     spotl0 = -spotl0;
%                   end
%               else
%                   if (rphase(t-1,:) >= pi/2 && rphase1 < pi/2) | (rphase(t-1,:) >= 3*pi/2 & rphase1 < 3*pi/2)
%                     spotl0 = -spotl0;
%                   end
%               end
%           end
          if Uint < 1. | abs(dxdt)<.1
              fprintf(1,'*** Uint < 0 or Ug->0, breaking: Uint=%d dxdt = %g (Xint,Yint)=(%g,%g)\n', Uint,dxdt,Xint,Yint);
              break
          end
                    %% make sure ray does not leave the domain where
          %% background fields are given:
          if Yint<lat(jmax) || Yint>lat(jmin)
            fprintf(1,'*** Yint>Ymax, breaking.  (Xint,Yint)=(%g,%g)\n'...
                    ,Xint,Yint);
            break
          end
             
          
%           if rem(t,day/(2*dt))==0
            fprintf(1,'t = %g spotk =  %g spotl = %g Ks = %g Ksint = %g  phase = %g \n ',t,spotk*rad*Cosint,spotl*rad,sqrt(Ks2int)*rad,sqrt(Ks2int)*rad,rphase1);
            fprintf(1,'t=%d (Xint,Yint)=(%g,%g) (spotk,spotl)=(%g,%g) \n'...
                    ,t, Xint,Yint,spotk*rad*Cosint,spotl*rad); 
%           end 
          
        end
        
        %writing out
         %the initial posistion and wavenumbers
        rphase0 = xm(i0)*real(spotk0(j0-jmin+1,i0)) + ym(j0)*real(spotl0(j0-jmin+1,i0)); % wave phase at t=0
        rphase10 = rphase1;

        if rphase10 > 2*pi 
              while rphase10 > 2*pi
                  rphase10 = rphase10 - 2*pi;
              end
        else if rphase10 < 0
              while rphase1 < 0
                  rphase10 = rphase10 + 2*pi;
              end
        end
        end
        rphase02pi = rphase10;

        
%         alL(1,:)=[0 0. 0. xm(i0) ym(j0) lon(i0) lat(j0) real(kk) real(spotl0(j0-jmin+1,i0))*rad imag(kk) imag(spotl0(j0-jmin+1,i0))*rad real(omega) imag(omega) rphase0 rphase02pi];
        alL(1,:)=[0 0. 0. xm(i0) ym(j0) lon(i0) lat(j0) real(kk) real(spotl0(j0-jmin+1,i0))*rad imag(kk) imag(spotl0(j0-jmin+1,i0))*rad rphase0 rphase02pi];
        t1 = t
        %for t = 1:Nsteps
        for t = 1:t1
           if rem(t*dt,day/kout)==0
               t
            if t<=size(timestep,1)
%                  alL(2*t*dt/day+1,:)=[timestep(t,:) locm(t,:) locgeo(t,:) rwnums(t,:) iwnums(t,:)  rsom(t,:)  isom(t,:) rphase(t,:) rphase2pi(t,:)];
                 alL(kout*t*dt/day+1,:)=[timestep(t,:) locm(t,:) locgeo(t,:) rwnums(t,:) iwnums(t,:) rphase(t,:) rphase2pi(t,:)];
            else
                break
            end
          end
        end
        ts = size(alL);
            if (uwnd_dim == 0);
                 fn_out = sprintf('../output/matlab/clim/rays/ray1d_critK2d_%s%d-%d_%dN_%dE_period%d_k%d_root%d',...
                              bgf,fyr1,lyr,lat0(ilat),lon0(ilon),period,kk,RR);
            end
            if (uwnd_dim == 1);
                fn_out = sprintf('../output/matlab/clim/rays/ray1d_realU_%s%d-%d_%dN_%dE_period%d_k%d_root%d',...
                             bgf,fyr1,lyr,lat0(ilat),lon0(ilon),period,kk,RR);
            end
        if ts(1)>1
            
            dlmwrite(fn_out, alL,'precision', '%.6f');
            sprintf( 'this is not skipped');
        else
            if exist(fn_out, 'file')
              % warningMessage = sprintf('Warning: file  exist:\n%s', fn_out);
              % uiwait(msgbox(warningMessage));

              delete (fn_out);
            end
        end 
            
        
        clear locm
        clear rwnums
        clear iwnums
        clear rw
        clear iw
        clear locgeo
        clear rsom
        clear isom
        clear rphase
        clear alL
       end
     end
    end
% end   %ilon
end

clearvars alL

% %------------------------------------------------
% %%%  !!! add time from the start of inegration (s) to the output !!!
% %fout = sprintf('../output/ray_loc%dE_%dN_period%d_k%d_root%d'...
% %                         ,lon0(),lat0(),period,kk,RR);
