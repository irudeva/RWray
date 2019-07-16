

close all;  clear all

addpath matlab matlab/spherepack
addpath /Applications/MATLAB_R2016a.app/m_map/

rad     = 6.371e6  ; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5 ; % rotation rate of Earth (rad/s)
dtr     = pi/180   ;
rtd     = 180/pi   ;

%% specify the number of points at the north/south pole to remove from analysis
j_pole=5;

% lat band
lat1 = [45,55,65 ];
lat2 = lat1+10;

%regions for temp
reg = ['C_Europe ';'E_Europe ';'W_Siberia';'E_Siberia';'W_China  ';'E_China  ';'W_America';'E_America' ];

lonreg1 = [ 5, 30, 70, 100,  80, 100, 235, 260 ];
lonreg2 = [25, 50, 90, 120, 100, 120, 255, 280 ];;
latreg1 = [40, 45, 50,  55,  25,  22,  45,  35 ];
latreg2 = [55, 60, 65,  70,  45,  42,  60,  50 ];


syr = 1979;
eyr = 2017;

ssn = ['DJF';'JJA';'MAM';'SON';'Jan';'Feb';'Mar';'Apr';'MAy';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec']

% fyr = 1979;
% lyr = 1985;

for yr = syr:eyr
    fyr = yr;
    lyr = fyr;

    % parameter setting: u wind
     % uwnd_dim = 0 to select 1 dim wind (zonally aver) 
     % uwnd_dim = 1 to select 2 dim wind (real wind)
    uwnd_dim = 1;

    % day = 24*60*60; % in s
    % dt  = 15*60;  % in s
    % integration_time=10*day;
    % Nsteps = round(integration_time/dt);
    % 
    % Periods ('Inf' for quasi-stationary waves) 
    %Periods=[ Inf 50 20 ]*day;
    % Periods=[ Inf ]*day;
    % 
    % freq=2*pi./Periods;
    % Nfr=length(freq);
    % 
    % kout = 24 % number of points per day in the output file

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
    %  lat0 = [65] ; %deg.N
    %  lon0 = [10] ; %deg E

    % smoothing before ray tracing MIGHT be a good idea...:
    % do_smooth_background_fields=1;

    % do_complex_tracing=0

    disp('-----------------------------------');


    %  disp(['k wavenumbers: ' num2str(k_wavenumbers(:)') '']) ;
    % disp(['Periods: ' num2str(Periods(:)'/day) ' days']) ;
    % disp(['Integration time: ', num2str(integration_time'/day),' days']) ;
    % disp(['dt: ' num2str(dt'/60) ' min']) ;
    % disp(['Stating point: ' num2str(lon0) ' E   ' num2str(lat0) ' N']) ;

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
    
%     ft     = sprintf('~/work/DATA/ERAint/Plev/erain.hgt_air_wind.monmean.%d.nc',iyr);
    ft     = sprintf('~/work/DATA/ERAint/Mslp_1deg/erain.mslp_air.monmean.%d.nc',iyr);
    disp(['Temp read from ' ft])
    ncid      = netcdf.open ( ft,'NC_NOWRITE' );
%     level     = netcdf.getVar (ncid,2);
%     ilev       = find(level==300);

    varid     = 4;
    [name,type,dimids,natts] = netcdf.inqVar(ncid,varid);
    if(name=='t2m'); 
        t_air = netcdf.getVar (ncid,varid,'short');
        scale = netcdf.getAtt(ncid,varid,'scale_factor');
        offset= netcdf.getAtt(ncid,varid,'add_offset');
        t_air = single(t_air)*scale+offset;   
%         t(:,:,:,iyr - syr+1) = squeeze(t_air(:,:,ilev,:));
        t(:,:,:,iyr - syr+1) = t_air(:,:,:);
         display(['size of t_air (',name,') = ',sprintf(' %d',size(t))]);
    else
        disp('Check t_air');
        exit
    end;



    % -------------------------------------------------------------------

    % Select time for the background field (climatology)

    TIMEbase = datenum(1900, 1, 1);
    date = datestr(double(time)/24 + TIMEbase); % where time is hours from Timebase

    formatdata = '01-%s-%d';
    %mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!

    for is = 5:size(ssn)
        bgf = ssn(is,:);
        disp(['season: ' bgf '']) ;
        if bgf =='DJF'
         clear mon;
         mon=[ 'Dec';'Jan';'Feb' ];  %%%!!!! for DEC -  year = year-1!!!
        elseif bgf=='JJA'
          clear mon;
          mon=[ 'Jun';'Jul';'Aug' ];  %%%!!!! for DEC -  year = year-1!!!
        elseif bgf=='MAM'
         clear mon;
         mon=[ 'Mar';'Apr';'May' ];  
        elseif bgf=='SON'
         clear mon;
         mon=[ 'Sep';'Oct';'Nov' ];  
        elseif bgf =='Jan'
            clear mon;
         mon=[ 'Jan'];  %%%!!!! for DEC -  year = year-1!!!
        elseif bgf=='Feb'
            clear mon;
         mon=[ 'Feb' ];  %%%!!!! for DEC -  year = year-1!!!
        elseif bgf=='Mar'
            clear mon;
         mon=[ 'Mar' ];  
        elseif bgf=='Apr'
            clear mon;
         mon=[ 'Apr' ];  
        elseif bgf =='May'
            clear mon;
         mon=[ 'May' ];  
        elseif bgf =='Jun'
            clear mon;
         mon=[ 'Jun'];  
        elseif bgf=='Jul'
            clear mon;
         mon=[ 'Jul' ]; 
        elseif bgf=='Aug'
            clear mon;
         mon=[ 'Aug' ];  
        elseif bgf=='Sep'
            clear mon;
         mon=[ 'Sep' ];  
        elseif bgf =='Oct'
            clear mon;
         mon=[ 'Oct' ]; 
        elseif bgf =='Nov'
            clear mon;
         mon=[ 'Nov' ]; 
        elseif bgf =='Dec'
            clear mon;
         mon=[ 'Dec' ]; 
        end

        in=1;
        fyr1 = fyr
        if bgf == 'DJF'
            fyr1 = fyr+1;
        end
        for yr = fyr1:lyr
            for imon = 1:size(mon(:,1))
              nyr = yr;
              if(mon(imon,:)=='Dec' & bgf=='DJF'); nyr=yr-1; end
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
              tssn(:,:,in) = zUbar(:,:,t1,t2);
             end
             end
         end
        end

        % Climatology
        u0=(squeeze(mean(ussn(:,:,:),3)))';
        [m,n]=size(u0);
        % v0=(squeeze(mean(vssn(:,:,:),3)))';
        % psi0=(squeeze(mean(psissn(:,:,:),3)))';
        t0=(squeeze(mean(tssn(:,:,:),3)))';

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
            t0=flipud(t0);
            disp('The data were flipped to be N->S')
        end


        % In Mercator projection

        xm=transpose(rad*lon*dtr);
        ym=rad*log((1+sin(dtr*lat))./cos(dtr*lat));
        ym(lat==90)=inf; % was using 1e10 with success here
        ym(lat==-90)=-inf;

        %  Smoothing the fields if do_smooth_background_fields=1
        % if do_smooth_background_fields
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
        % end


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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %masking
        ind1 = find(UbarM<3);
        ind2 = find(BetaM<0);
        %ind3 = find(isnan(BetaM));

        BetaM1 = BetaM;
        BetaM1(ind1) = NaN;
        BetaM1(ind2) = NaN;
        UbarM1 = UbarM;
        UbarM1(ind1) = NaN;
        UbarM1(ind2) = NaN;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ks 

        % sKs2=abs(rdivide(sBetaM1,sUbarM1-rdivide(omega,spotk0)));
        % stationary waves:
        Ks=real (sqrt(rdivide(BetaM1,UbarM1)))*rad;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotting Ks maps

%         clf;
%         axes('position',[.1 .2 .8 .7]);
%         m_proj('mercator','lon',[0 360],'lat',[-85 85]);
% 
%         % [cs,h]= m_contourf(lon,lat(6:176),(real(sqrt(rdivide(sBetaM1,sUbarM1)))*rad));
%         [cs,h]= m_contourf(lon,lat,Ks,[0:1:20]);
% 
%         % hold on;
% 
%          m_coast('linewidth',1,'color','r');;
%          m_grid;
% 
% %          title(sprintf('Ks, %s%d-%d',bgf,fyr1,lyr),'fontsize',14);
%          title(sprintf('Ks, %s%d',bgf,fyr1),'fontsize',14);
%          ax=m_contfbar([.2 .8],-.1,cs,h);
%          set(ax,'fontsize',12)
% 
% %         figout     = sprintf('../output/Ks/Ks.%s%d-%d.png',bgf,fyr1,lyr); 
%         figout     = sprintf('../output/Ks/Ks.%s%d.png',bgf,fyr1); 
%         saveas(gcf,figout)


       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time series
        
        for iband = 1:3
            
            ilat1=find(lat==lat1(iband));
            ilat2=find(lat==lat2(iband));
            Kwhite(iband,is,fyr1-syr+1) = sum(sum(isnan(Ks(ilat2:ilat1,:))));
            Ktotal(iband,is,fyr1-syr+1) = numel(Ks(ilat2:ilat1,:));
      
        end
        
        for ireg = 1:length(lonreg1)
            
            ilon1=find(lon==lonreg1(ireg));
            ilon2=find(lon==lonreg2(ireg));
            ilat1=find(lat==latreg1(ireg));
            ilat2=find(lat==latreg2(ireg));
            treg(ireg,is,fyr1-syr+1) = mean2(t0(ilat2:ilat1,ilon1:ilon2));
      
        end
 
        
    end   %bgf
end  %yr



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting time series
% 
% f = figure;
% 
% for is = 5:size(ssn)
%     clf(f);
%     bgf = ssn(is,:);
%     disp(['season  plot: ' bgf '']) ;
%     pK = plot([syr:eyr],100*squeeze(Kwhite(:,is,:)/ Ktotal(1,is,fyr1-syr+1)));
%     ylim([0 80])
%     legend(sprintf('%d-%d',lat1(1),lat2(1)),sprintf('%d-%d',lat1(2),lat2(2)),sprintf('%d-%d',lat1(3),lat2(3)));
%     title(sprintf('Ks white, %%, %s',bgf),'fontsize',14);
% %     hold on;
%     yyaxis right
%     pT = plot([syr:eyr],squeeze(treg(:,is,:)),'--');
%     pT(1).Color = 'r';
%     pT(2).Color = 'b';
%     pT(3).Color = [0.9290,0.6940,0.1250];
%     pT(4).Color = 'g';
%     ylabel('Temperature');
%     ylim([-40 30])
% 
%     lgd = legend(sprintf('%d-%d',lat1(1),lat2(1)),sprintf('%d-%d',lat1(2),lat2(2)),sprintf('%d-%d',lat1(3),lat2(3)),reg(1,:),reg(2,:),reg(3,:),reg(4,:),'Location','northwest');
% %     lgd.NumColumns = 2
%     lgd.Position = [0.7482 0.7012 0.1313 0.2000];
%     
% %  add table
% %     f = figure;
%     t = uitable(f);
% 
%     tabl(1,:) ={'region',reg(1,:),reg(2,:),reg(3,:),reg(4,:)};
%     for ir = 1:length(lat1(1,:))
%         tabl(ir+1,1) = {sprintf('%d-%d',lat1(ir),lat2(ir))};
%         for ic = 1:length(reg(:,1))
%           corr = corrcoef(100*squeeze(Kwhite(ir,is,:)/ Ktotal(ir,is,fyr1-syr+1)),squeeze(treg(ic,is,:)))
%           tabl(ir+1,ic+1) = {corr(1,2)};
%         end
%     end
% 
%     t.Data = tabl;
%     t.Position = [20 20 410 105];
% 
%     figout     = sprintf('../output/Ks/Ks_Treg_corr.%s%d-%d.png',bgf,syr,eyr); 
%     saveas(gcf,figout);
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  table

for ir = 1:length(lat1(1,:))
        region(ir,:) = {sprintf('%d-%d',lat1(ir),lat2(ir))};
        
        for is=5:16
        
        for ic = 1:length(reg(:,1))
            corr = corrcoef(100*squeeze(Kwhite(ir,is,:)/ Ktotal(ir,is,fyr1-syr+1)),squeeze(treg(ic,is,:)));
%           tbl(ir+1,ic+1) = {corr(1,2)};
            corrtbl(is,ir,ic) = corr(1,2);
        end
        end
        
end

varNames = {'region   '};
    for ic = 1:length(reg(:,1))
        varNames(1,ic+1) = {reg(ic,:)};
    end

 
 for is=5:16
     varNames(1,1) = {ssn(is,:)};
     tmp =  squeeze(corrtbl(is,:,:));
     T = table(region, tmp(:,1));
     for ir =2:length(reg(:,1))
         T(:,ir+1)= table(tmp(:,ir));
     end
     T.Properties.VariableNames = varNames;
     
     writetable(T,'../output/Ks/Ks_Treg_corr.xls','Sheet',1,'Range',sprintf('A%d:I%d',(is-5)*5+1,(is-5)*5+5))
 end

