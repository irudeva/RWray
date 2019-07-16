

close all;  clear all

addpath matlab matlab/spherepack
addpath /Applications/MATLAB_R2016a.app/m_map/

rad     = 6.371e6  ; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5 ; % rotation rate of Earth (rad/s)
dtr     = pi/180   ;
rtd     = 180/pi   ;

%% specify the number of points at the north/south pole to remove from analysis
j_pole=5;


syr = 1979;
eyr = 1979;

%ssn = ['DJF';'JJA';'MAM';'SON';'Jan';'Feb';'Mar';'Apr';'MAy';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec']

% parameter setting: u wind
 % uwnd_dim = 0 to select 1 dim wind (zonally aver)
 % uwnd_dim = 1 to select 2 dim wind (real wind)
uwnd_dim = 1;


disp('-----------------------------------');


%  disp(['k wavenumbers: ' num2str(k_wavenumbers(:)') '']) ;
% disp(['Periods: ' num2str(Periods(:)'/day) ' days']) ;
% disp(['Integration time: ', num2str(integration_time'/day),' days']) ;
% disp(['dt: ' num2str(dt'/60) ' min']) ;
% disp(['Stating point: ' num2str(lon0) ' E   ' num2str(lat0) ' N']) ;

% disp('-----------------------------------');
disp('');
% Read data

for iyr = syr:eyr

    fin     = sprintf('~/work/DATA/ERAint/Plev/erain.hgt_air_wind.6h.%d.nc',iyr);
    disp(['Data read from ' fin])
    ncid      = netcdf.open ( fin,'NC_NOWRITE' );

    [name,type,dimids,natts] = netcdf.inqVar(ncid,2);
    if(name == 'level');
        lev     = netcdf.getVar (ncid,2);
    %         levin   = lev;
    %         nlev    = size(lev,1);
    %          display(['size of level = ',sprintf(' %d',nlev)]);
    %    ilev  = find(lev==300);

    else
        disp('Check lev');
        exit
    end;



    varid     = 7;
    [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name=='u');
        uwnd     = netcdf.getVar (ncid,varid,'short');
        scale = netcdf.getAtt(ncid,varid,'scale_factor');
        offset= netcdf.getAtt(ncid,varid,'add_offset');
        uwnd = single(uwnd)*scale+offset;
%         u(:,:,:) = squeeze(uwnd(:,:,ilev,:));
        u = uwnd;
        u_1MS = dset['v'].resample(time='1MS').mean('time')

         display(['size of uwnd (',name,') = ',sprintf(' %d',size(u))]);
    else
        disp('Check u');
        exit
    end;
    varid     = 8;
    [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name=='v');
        vwnd     = netcdf.getVar (ncid,varid,'short');
        scale = netcdf.getAtt(ncid,varid,'scale_factor');
        offset= netcdf.getAtt(ncid,varid,'add_offset');
        vwnd = single(vwnd)*scale+offset;
%         v(:,:,:) = squeeze(vwnd(:,:,ilev,:));
        v = uwnd;
         display(['size of vwnd (',name,') = ',sprintf(' %d',size(v))]);
    else
        disp('Check v');
        exit
    end;
%     varid     = 5;
%     [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
%      if(name=='sf');
%         %disp('reading streamfunction');
%         psi1     = netcdf.getVar (ncid,varid);
%         scale = netcdf.getAtt(ncid,varid,'scale_factor');
%         offset= netcdf.getAtt(ncid,varid,'add_offset');
%         psi(:,:,:,iyr - fyr+1) = single(psi1)*scale+offset;
%          display(['size of streamfunction (',name,') = ',sprintf(' %d',size(psi))]);
%     else
%         disp('Check sf');
%         exit
%     end;

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

    [name,type,dimids,natts] = netcdf.inqVar(ncid,3);
    if(name == 'time');
        time = netcdf.getVar (ncid,3);
         display(['size of time = ',sprintf(' %d',size(time))]);
    else
        disp('Check time');
        exit
    end;
%     netcdf.close(ncid);

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
        zVbar = v;
    end

    % -------------------------------------------------------------------

    % Select time for the background field (climatology)

    TIMEbase = datenum(1900, 1, 1);
    date = datestr(double(time)/24 + TIMEbase); % where time is hours from Timebase

    formatdata = '%d-%s-%d %d:00:00';
    %mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!

   for it =1:size(time,1)
       disp(['Time: ' num2str(it) ]) ;
       for il = 1:size(lev)
%         disp(['Level: ' num2str(lev(il)) ]) ;
        ussn = zUbar(:,:,il,it);
        vssn = zVbar(:,:,il,it);


        % Climatology
        % u0=(squeeze(mean(ussn(:,:,:),3)))';
        u0=(squeeze(ussn))';
        [m,n]=size(u0);
        v0=(squeeze(vssn))';

        % Check that data are from N to S
        if lat(1)<lat(end)
            if size(u0,1)~=length(lat) || size(u0,2)~=length(lon)
              disp('Flipud: u0 array mismatched to dimensions')
              exit
            end
            lat=flipud(lat);
            u0=flipud(u0);
            % v0=flipud(v0);
            disp('The data were flipped to be N->S')
        end


        % In Mercator projection

        xm=transpose(rad*lon*dtr);
        ym=rad*log((1+sin(dtr*lat))./cos(dtr*lat));
        ym(lat==90)=inf; % was using 1e10 with success here
        ym(lat==-90)=-inf;

        %  Smoothing the fields if do_smooth_background_fields=1
        % if do_smooth_background_fields
        u0=zfltr(u0,1,10,1);
        v0=zfltr(v0,1,10,1);
         % end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%  Solving for UbarM, VbarM, and even psiM

        [X,Y]=meshgrid(lon,lat);
        [Lon,Lat]=meshgrid(lon*dtr,lat*dtr);

        UbarM=u0./cos(Lat);
        % VbarM=v0./cos(Lat);
        wind = sqrt(u0.^2+v0.^2);

        UbarM(Y==90)=inf; % was using 1e10 with success here
        UbarM(Y==-90)=-inf;
        % VbarM(Y==90)=inf; % was using 1e10 with success here
        % VbarM(Y==-90)=-inf;
        % PsiM(Y==90)=inf; % was using 1e10 with success here
        % PsiM(Y==-90)=-inf;

        if (uwnd_dim == 1);
            UbarM_2d = UbarM;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %                                                 continue
        end
        cdy(1:2,:) = nan ;
        cdy(nlat-1:nlat,:) = nan ;
        cdy=cdy.*d;
        for j=4:nlat-3
          cdy2(j,:)=(mean(cdy(j-1:j,:))-mean(cdy(j:j+1,:)))/(mean(ym(j- ...
                                                            1:j))-mean(ym(j:j+1)));
    %       continue
        end
        cdy2(1:3,:) = nan ;
        cdy2(nlat-2:nlat,:) = nan ;

        BetaM=trm1-cdy2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %masking
    %     ind1 = find(UbarM<3);
        ind1 = find(UbarM<5);
        ind2 = find(BetaM<0);
        %ind3 = find(isnan(BetaM));

        BetaM1 = BetaM;
        BetaM1(ind1) = NaN;
        BetaM1(ind2) = NaN;
        UbarM1 = UbarM;
        UbarM1(ind1) = NaN;
        UbarM1(ind2) = NaN;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Ks

        Ks(:,:,il,it)=(real(sqrt(rdivide(BetaM1,UbarM1)))*rad)';
        BetaM_3d(:,:,il,it) = BetaM1';

      end

   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing Ks to netcdf

    fout = sprintf('../output/Ks/erain.Ks_wind.6h.%d.nc',iyr)
    ncoutid = netcdf.create(fout,'NETCDF4');

    % define dimensions
    dimid_lon = netcdf.defDim(ncoutid,'longitude',nlon);
    dimid_lat = netcdf.defDim(ncoutid,'latitude',nlat);
    dimid_lev = netcdf.defDim(ncoutid,'level',size(lev,1));
    dimid_t = netcdf.defDim(ncoutid,'time',netcdf.getConstant('NC_UNLIMITED'));

    % define variables
    varid_lon = netcdf.defVar(ncoutid,'longitude','NC_FLOAT',dimid_lon);
    varid_lat = netcdf.defVar(ncoutid,'latitude','NC_FLOAT',dimid_lat);
    varid_lev = netcdf.defVar(ncoutid,'level','NC_INT',dimid_lev);
    varid_t = netcdf.defVar(ncoutid,'time','NC_INT',dimid_t);
%     varid_u = netcdf.defVar(ncoutid,'u','NC_SHORT',[dimid_lon dimid_lat dimid_lev dimid_t]);
%     varid_v = netcdf.defVar(ncoutid,'v','NC_SHORT',[dimid_lon dimid_lat dimid_lev dimid_t]);
    varid_Ks = netcdf.defVar(ncoutid,'Ks','NC_FLOAT',[dimid_lon dimid_lat dimid_lev dimid_t]);
    varid_BetaM = netcdf.defVar(ncoutid,'BetaM','NC_FLOAT',[dimid_lon dimid_lat dimid_lev dimid_t]);

    % attribute
    varid     = 0;
   [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name == 'longitude');
        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_lon);
        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_lon);
    end

    varid     = 1;
   [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name == 'latitude');
        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_lat);
        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_lat);
    end

        varid     = 2;
   [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name == 'level');
        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_lev);
        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_lev);
    end

    varid     = 3;
   [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name == 'time');
        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_t);
        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_t);
        netcdf.copyAtt(ncid,varid,'calendar',ncoutid,varid_t);
    end

    netcdf.putAtt(ncoutid,varid_Ks,'long_name','total stationary wavenumber');
    netcdf.putAtt(ncoutid,varid_Ks,'equation','Ks*Rad');
    Ks(isnan(Ks)) = -32767.;
    netcdf.defVarFill(ncoutid,varid_Ks,false,-32767.);

    netcdf.putAtt(ncoutid,varid_BetaM,'long_name','meridional gradient of absolute vorticity ');
    netcdf.putAtt(ncoutid,varid_BetaM,'equation','m**-1 s**-1');
    netcdf.defVarFill(ncoutid,varid_BetaM,false,-32767.);
    scale = 10^-12;
    netcdf.putAtt(ncoutid,varid_BetaM,'scale_factor',scale);
    BetaM_3d = BetaM_3d/scale;
    BetaM_3d(isnan(BetaM_3d)) = -32767.;


%      varid     = 7;
%     [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
%     if(name=='u');
%        netcdf.copyAtt(ncid,varid,'scale_factor',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'add_offset',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'_FillValue',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'standard_name',ncoutid,varid_u);
%     end if
%
%      varid     = 8;
%     [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
%     if(name=='v');
%        netcdf.copyAtt(ncid,varid,'scale_factor',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'add_offset',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'_FillValue',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'units',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'long_name',ncoutid,varid_u);
%        netcdf.copyAtt(ncid,varid,'standard_name',ncoutid,varid_u);
%     end if

    % leave define mode and enter data mode
    netcdf.endDef(ncoutid);

    % write data
    netcdf.putVar(ncoutid,varid_t,0,length(time),time);
%     netcdf.putVar(ncoutid,varid_t,0,5,time(1:5));
    netcdf.putVar(ncoutid,varid_lev,lev);
    netcdf.putVar(ncoutid,varid_lat,lat);
    netcdf.putVar(ncoutid,varid_lon,lon);
    netcdf.putVar(ncoutid,varid_Ks,Ks);
    netcdf.putVar(ncoutid,varid_BetaM,BetaM_3d);

    % close output and input netcdf
    netcdf.close(ncoutid);
    netcdf.close(ncid);

    clear Ks
    clear BetaM_3d


end


%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % plotting Ks maps
%
%         clf;
%         axes('position',[.1 .2 .8 .7]);
% %         m_proj('mercator','lon',[0 360],'lat',[-85 85]);
%         m_proj('stereographic','lat',90,'long',0,'radius',75);
%         m_grid('xtick',12,'tickdir','out','ytick',[30 40 50 60 70 80],'linest','-');
%
%         % [cs,h]= m_contourf(lon,lat(6:176),(real(sqrt(rdivide(sBetaM1,sUbarM1)))*rad));
% %         [cs,h]= m_contourf(lon,lat,wind,[-100:10:200]);
% %         [cs,h]= m_contourf(lon,lat,MCI,[-5:.2:5]);
% %         [cs,h]= m_contourf(lon,lat,wind,[0:5:70],'edgecolor','none');
%          Ks(:,361) = Ks(:,1);
%          lon361 = lon;
%          lon361(361) = 360.;
%         [cs,h]= m_contourf(lon361,lat,Ks,[0:4:20],'edgecolor','none');
%          m_coast('linewidth',1,'color','r');
%          m_grid;
% %          title(sprintf('Ks, %s%d, days: %d-%d',bgf,iyr,date1(id),date2(id)),'fontsize',14);
% %          title(sprintf('Ks, %d%s%d',iyr,bgf,date1(id)),'fontsize',14);
%          title({sprintf('Wind speed and Ks, 00:00 UTC %d%s%d',date1(id),bgf,iyr),''},'fontsize',14);
%          ax=m_contfbar([.2 .8],-.1,cs,h);
%          title(ax,'Ks');
%          set(ax,'fontsize',12)
%         hold on;
%
% %         [cs1,h1]= m_contour(lon,lat,Ks,[0:4:20],'edgecolor','k');
%         [cs1,h1]= m_contour(lon,lat,wind,[20:5:70],'edgecolor','k');
%         h1.LineWidth = 1.
%
%
% %         figout     = sprintf('../output/Ks/Ks.%s%d.%d_%d.png',bgf,iyr,date1(id),date2(id));
%         figout     = sprintf('../output/Ks/Ks.%d%s%d.png',iyr,bgf,date1(id));
% %         figout     = sprintf('../output/Ks/MCI23.%d%s%d.png',iyr,bgf,date1(id));
%         saveas(gcf,figout)


       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time series

%         for iband = 1:3
%
%             ilat1=find(lat==lat1(iband));
%             ilat2=find(lat==lat2(iband));
%             Kwhite(iband,is,yr-syr+1) = sum(sum(isnan(Ks(ilat2:ilat1,:))));
%             Ktotal(iband,is,yr-syr+1) = numel(Ks(ilat2:ilat1,:));
%
%         end
%
%         for ireg = 1:length(lonreg1)
%
%             ilon1=find(lon==lonreg1(ireg));
%             ilon2=find(lon==lonreg2(ireg));
%             ilat1=find(lat==latreg1(ireg));
%             ilat2=find(lat==latreg2(ireg));
%             treg(ireg,is,fyr1-syr+1) = mean2(t0(ilat2:ilat1,ilon1:ilon2));
%
%         end


%     end ; % days
% end  %yr
