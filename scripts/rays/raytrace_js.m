%%%%%%%%%%%%%%%%%%%%%%%%
%%  Script to do 2-D smoothed field ray trace for the full field
%%  for El Nino JFM climatology only
%%  There are 3 background states to be used, possibly
%%  Everything is done with units this time
%%  All fields needed are solved for in the grid and those grids are
%%  then interpolated in the actual ray tracing.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  The following are needed (based on Karoly, 1983):
%%  1)  BetaM_js--this can be solved for a number of ways, as Karoly did 
%%  (and also Hoskins and Karoly, 1981), where BetaM_js=2*Omega*cos2(Lat)/r
%%         - d/dy[(1/cos2(Lat))*d/dy(cos2(Lat)*UbarM_js)]
%%     or as Hoskins and Ambrizzi specify (2 ways: w/ and w/o Mercator) 
%%
%%  2)  UbarM_js--per Karoly, this is u/cos(Lat)
%%  3)  VbarM_js--per Karoly, this is v/cos(Lat)
%%  4)  x--r*Lat, where Lat=degrees*pi/180
%%  5)  y--r*log[(1+sin(lat))/cos(lat)]
%%  6)  qbar_js--this is (d2/dx+d2/dy)psi/cos2(lat) + 2*Omega*sin(lat)
%%  7)  dqbar_js/dy--this is BetaM_js
%%  8)  dqbar_js/dx--this is -d/dx[(d2/dx+d2/dy)psi/cos2(lat)]
%%  9)  d2qbar_js/dxx
%%  10) d2qbar_js/dxy
%%  11) d2qbar_js/dyy
%%  12) dUbarM_js/dx
%%  13) dUbarM_js/dy
%%  14) dVbarM_js/dx
%%  15) dVbarM_js/dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  From these we can interpolate the fields and calculate values needed
%%  for ray tracing, namely:
%%  Ks--this is (abs(BetaM_js/Um))^0.5 for zonally symmetric flows
%%     Kw, those with non-zero frequency are (abs(BetaM_js/(Um-w/k))^0.5
%%     For 2-D flow, Ks might be (abs(BetaM_js/Um + dqbar_js/dx/Vm))^0.5 or some such
%%  The other way is to calculate it for stationary waves straight from the
%%    wave numbers
%%  A) Ks=sqrt(k^2+l^2)
%%  B) dk/dt=-k*dUbarM_js/dx -l*dVbarM_js/dx +(d2qbar_js/dxy*k-d2qbar_js/dxx*l)/Ks^2
%%  C) dl/dt=-k*dUbarM_js/dy -l*dVbarM_js/dy +(d2qbar_js/dyy*k-d2qbar_js/dxy*l)/Ks^2
%%  D) dx/dt=ug=UbarM_js+{(k^2-l^2)*dqbar_js/dy - 2*k*l*dqbar_js/dx}/Ks^4
%%  E) dy/dt=vg=VbarM_js+{2*k*l*dqbar_js/dy - (k^2-l^2)*dqbar_js/dx}/Ks^4

%%  Start loading streamfunction, and winds

% $$$ loaddap('http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Intrinsic/.PressureLevel/.psi/T/%28Jan-Mar%29VALUES/P/%28200%29VALUES%5BT%5Daverage/dods')
% $$$ psi=squeeze(psi.psi);
% $$$ 
% $$$ loaddap('http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Intrinsic/.PressureLevel/.u/T/%28Jan-Mar%29VALUES/P/%28200%29VALUES%5BT%5Daverage/dods')
% $$$ u=squeeze(u.u);
% $$$ 
% $$$ loaddap('http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Intrinsic/.PressureLevel/.v/T/%28Jan-Mar%29VALUES/P/%28200%29VALUES%5BT%5Daverage/dods')
% $$$ v=squeeze(v.v);

ncid = netcdf.open('../output/zon_var/test.u.nc');
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);
for i=0:numvars-1
[name,xtype,dimids,natts] = netcdf.inqVar(ncid,i);
u=netcdf.getVar(ncid,i);u=u';
end
Y=netcdf.getVar(ncid,1);
X=netcdf.getVar(ncid,0);

ncid = netcdf.open('../output/zon_var/test.v.nc');
v=netcdf.getVar(ncid,2);v=v';
ncid = netcdf.open('../output/zon_var/test.psi.nc');
psi=netcdf.getVar(ncid,2);
  scale = netcdf.getAtt(ncid,2,'scale_factor');
  offset= netcdf.getAtt(ncid,2,'add_offset');
  psi = single(psi)*scale+offset;
psi=psi';

x=transpose(X); 
y=Y;  
Y=y*ones(1,480);
X=ones(241,1)*x;
lat=Y*pi/180;
Lat(:,:,1)=lat;

f=(2*7.2925e-5)*sin(Lat);
r=6.37e6;

%%%%%%%%%%%%%%%%%%%%
%%  Getting the x and y Mercator coordinates

xx=r*x*pi/180;
yy=r*double(log((1+sin(y*pi/180))./cos(y*pi/180)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Smoothing the fields

for SM=2
if SM==2
u200=zfltr(u,1,10,1);
v200=zfltr(v,1,10,1);
psi200=zfltr(psi,1,10,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for UbarM_js, VbarM_js, and even PsiM_js

%%%  Ks= {[2*omega - ((1/cos(lat))(d/d(lat))^2)*((1+cos2*Lat)/2)*U]/U}^0.5

UbarM_js=u200./cos(Lat);
VbarM_js=v200./cos(Lat);
PsiM_js=psi200./cos(Lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for BetaM_js; NOTE that cos2(Lat)=(1+cos(2*Lat))/2

a=(2*7.2925e-5);
b=(1+cos(2*Lat))/2;
trm1=a*b/r;
c=b.*UbarM_js;
d=1./b;
cdy=c;cdy2=c;
for i=2:240
     cdy(i,:)=(mean(c(i-1:i,:))-mean(c(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
cdy=cdy.*d;
for i=3:239
     cdy2(i,:)=(mean(cdy(i-1:i,:))-mean(cdy(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end

BetaM_js=trm1-cdy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for qbar_js
%%%%  The larger structure is to make sure have entire zonal field

tempPSI=[PsiM_js(:,477:480) PsiM_js PsiM_js(:,1:4)];
tempxx=[xx(477:480) xx xx(1:4)];

px1=tempPSI;px2=px1;py1=px1;py2=px1;
for i=2:240
     py1(i,:)=(mean(tempPSI(i-1:i,:))-mean(tempPSI(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
for i=3:239
     py2(i,:)=(mean(py1(i-1:i,:))-mean(py1(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
for i=2:487
     px1(:,i)=(mean(tempPSI(:,i-1:i)')'-mean(tempPSI(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end
for i=3:486
     px2(:,i)=(mean(px1(:,i-1:i)')'-mean(px1(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end

dd=[d(:,477:480) d d(:,1:4)];
tempf=[f(:,477:480) f f(:,1:4)];

tempqbar_js=tempf+dd.*(px2+py2);
qbar_js=tempqbar_js(:,5:484);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dqbar_js/dx, and d2qbar_js/dx2

px1=NaN*ones(241,488); px2=px1;
for i=4:485
     px1(:,i)=(mean(tempqbar_js(:,i-1:i)')'-mean(tempqbar_js(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end
for i=5:484
     px2(:,i)=(mean(px1(:,i-1:i)')'-mean(px1(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end

dqbar_jsdx=px1(:,5:484);
d2qbar_jsdx2=px2(:,5:484);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dqbar_js/dy (should be same as BetaM_js), and d2qbar_js/dy2

py1=NaN*ones(241,488); py2=py1;
for i=4:238
     py1(i,:)=(mean(tempqbar_js(i-1:i,:))-mean(tempqbar_js(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
for i=5:237
     py2(i,:)=(mean(py1(i-1:i,:))-mean(py1(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end

dqbar_jsdy=py1(:,5:484);
d2qbar_jsdy2a=py2(:,5:484);

%%%%  Alternately for d2qbar_js/dy2

py2=NaN*ones(241,480);
for i=4:238
     py2(i,:)=(mean(BetaM_js(i-1:i,:))-mean(BetaM_js(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end

d2qbar_jsdy2=py2;

sum(sum(d2qbar_jsdy2(6:236,6:474)-d2qbar_jsdy2a(6:236,6:474)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for d2qbar_js/dxy (once);  COULD ALSO DO dx of either BetaM_js or
%%%%  dqbar_jsdy

py1=NaN*ones(241,480);
for i=2:240
     py1(i,:)=(mean(dqbar_jsdx(i-1:i,:))-mean(dqbar_jsdx(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end

d2qbar_jsdxy=py1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dUbarM_js/dx and dUbarM_js/dy

tempU=[UbarM_js(:,480) UbarM_js UbarM_js(:,1)];

px1=NaN*ones(241,482);py1=px1;
for i=2:240
     py1(i,:)=(mean(tempU(i-1:i,:))-mean(tempU(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
for i=2:481
     px1(:,i)=(mean(tempU(:,i-1:i)')'-mean(tempU(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end

dUbarM_jsdx=px1(:,2:481);
dUbarM_jsdy=py1(:,2:481);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for dVbarM_js/dx and dVbarM_js/dy

tempV=[VbarM_js(:,480) VbarM_js VbarM_js(:,1)];

px1=NaN*ones(241,482);py1=px1;
for i=2:240
     py1(i,:)=(mean(tempV(i-1:i,:))-mean(tempV(i:i+1,:)))/(mean(yy(i-1:i))-mean(yy(i:i+1)));
end
for i=2:481
     px1(:,i)=(mean(tempV(:,i-1:i)')'-mean(tempV(:,i:i+1)')')/(mean(tempxx(i-1:i))-mean(tempxx(i:i+1)));
end

dVbarM_jsdx=px1(:,2:481);
dVbarM_jsdy=py1(:,2:481);

%figure(111)
%set(gcf,'PaperPosition',[0.5 0.5 7.5 10])
%subplot(2,1,1)
%m_proj('equidistant','lon',[0 358],'lat',[-88 88]);
%m_pcolor(X,Y,u);shading flat
%colorbar
%m_coast('linewidth',1,'color','k');
%m_grid
%title('\psi')
%subplot(2,1,2)
%m_proj('equidistant','lon',[0 358],'lat',[-88 88]);
%m_pcolor(X,Y,u200);shading flat
%m_coast('linewidth',1,'color','k');
%m_grid
%colorbar
%title('Smoothed \psi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving the system for different forcing sites

%%%% Irina Rudeva wants 30N 30E
xx=double(xx);
frcy=[77 77 77 70 70 70 84 84 84];  %  On reduced grid (5:237) 30N is 77
                            %  (trying 35.25 and 24.75, too)
frcx=[41 34 48 41 34 48 41 34 48];
FRR=[3];
jmin=5;jmax=237;
run_record=[];

for fr=1:1  %% IR
     frx=frcx(fr);
     fry=frcy(fr);

[subxx,subyy]=meshgrid(xx,yy(5:237));
[subyy_interp2,subxx_interp2]=meshgrid(yy(5:237),xx);
[XX,YY]=meshgrid(x,y(5:237));
sBetaM_js=BetaM_js(5:237,:);
sUbarM_js=UbarM_js(5:237,:);
sVbarM_js=VbarM_js(5:237,:);
sqbar_js=qbar_js(5:237,:);
sdqbar_jsdy=dqbar_jsdy(5:237,:);
sdqbar_jsdx=dqbar_jsdx(5:237,:);
sd2qbar_jsdx2=d2qbar_jsdx2(5:237,:);
sd2qbar_jsdy2=d2qbar_jsdy2(5:237,:);
sd2qbar_jsdxy=d2qbar_jsdxy(5:237,:);
sdUbarM_jsdx=dUbarM_jsdx(5:237,:);
sdUbarM_jsdy=dUbarM_jsdy(5:237,:);
sdVbarM_jsdx=dVbarM_jsdx(5:237,:);
sdVbarM_jsdy=dVbarM_jsdy(5:237,:);

%%%%%%%%%%%%%
%%%%  Estimating the initial Ks from the forcing site for a specific
%%%%  BetaM_js and UbarM_js


rrr=[1];
for i=1:length(rrr)
dt=3600/rrr(i);
for kkz=1:length(FRR)
  kk=FRR(kkz); 
for RR=1:3
for hh=1
  spotk=kk/r/cos(Lat(fry+4,frx));
  subk=[fr kk RR]
  ytt=yy(5:237);yi=ytt(fry);xi=xx(frx);
%  Ks=sqrt(abs(sBetaM_js(fry,frx)/sUbarM_js(fry,frx)));
%  spotl=real((Ks^2-spotk^2)^0.5);
Bint=sBetaM_js(fry,frx);Uint=sUbarM_js(fry,frx);Vint=sVbarM_js(fry,frx);
qxint=sdqbar_jsdx(fry,frx);qyint=sdqbar_jsdy(fry,frx);
cz(1)=Vint;
cz(2)=Uint*spotk;
cz(3)=Vint*spotk^2+qxint;
cz(4)=Uint*spotk^3-qyint*spotk;
tl=roots(cz);
     spotl=tl(RR);
     Ks=(spotl^2+spotk^2)^0.5;

%%%%%%%%%%%%%%%%
%%  Breaking the run if the starting l wavenumber is too big

if abs(spotl*r*cos(Lat(fry+4,1)))>30
  abs(spotl*r*cos(Lat(fry+4,1)))
break
elseif imag(spotl)~=0
break
elseif spotl<0
end

%%%%%%%%%%%
%%  Starting the loop with the above initial k,l, and Ks

  for t=1:360*rrr(i)
   if rem(t,60)==0
    t 
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Interpolating the fields to the current spot

Uint=interp2(subyy_interp2,subxx_interp2,sUbarM_js',yi,xi,'spline');
Bint=interp2(subyy_interp2,subxx_interp2,sBetaM_js',yi,xi,'spline');
Vint=interp2(subyy_interp2,subxx_interp2,sVbarM_js',yi,xi,'spline');
qint=interp2(subyy_interp2,subxx_interp2,sqbar_js',yi,xi,'spline');
qyint=interp2(subyy_interp2,subxx_interp2,sdqbar_jsdy',yi,xi,'spline');
qxint=interp2(subyy_interp2,subxx_interp2,sdqbar_jsdx',yi,xi,'spline');
qyyint=interp2(subyy_interp2,subxx_interp2,sd2qbar_jsdy2',yi,xi,'spline');
qxxint=interp2(subyy_interp2,subxx_interp2,sd2qbar_jsdx2',yi,xi,'spline');
qxyint=interp2(subyy_interp2,subxx_interp2,sd2qbar_jsdxy',yi,xi,'spline');
Uyint=interp2(subyy_interp2,subxx_interp2,sdUbarM_jsdy',yi,xi,'spline');
Uxint=interp2(subyy_interp2,subxx_interp2,sdUbarM_jsdx',yi,xi,'spline');
Vyint=interp2(subyy_interp2,subxx_interp2,sdVbarM_jsdy',yi,xi,'spline');
Vxint=interp2(subyy_interp2,subxx_interp2,sdVbarM_jsdx',yi,xi,'spline');
Cosint=interp2(subyy_interp2,subxx_interp2,cos(Lat(jmin:jmax,:))',yi,xi,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Solving for the changes

dkdt=-spotk*Uxint-spotl*Vxint+(qxyint*spotk-qxxint*spotl)/Ks^2;
dldt=-spotk*Uyint-spotl*Vyint+(qyyint*spotk-qxyint*spotl)/Ks^2;
dxdt=Uint+((spotk^2-spotl^2)*qyint-2*spotk*spotl*qxint)/Ks^4;
dydt=Vint+(2*spotk*spotl*qyint-(spotk^2-spotl^2)*qxint)/Ks^4;

%%%%%%%%%%%%%%
%%  Updating the changes

xi=xi+dxdt*dt;
if xi>=max(max(subxx))
     xi=xi-max(max(subxx));
end
yi=yi+dydt*dt;
spotl=spotl+dldt*dt;
spotk=spotk+dkdt*dt;
Ks=(spotk^2+spotl^2)^0.5;

%%%%%%%%%%%%%%%
%%  Finding the location

Yint=interp2(subyy_interp2,subxx_interp2,YY',yi,xi,'spline');
Xint=interp2(subyy_interp2,subxx_interp2,XX',yi,xi,'spline');

%%%%%%%%%%%%%%
%%  Storing

adj=r*Cosint;
trl(t,:)=[xi yi];
nums(t,:)=[spotl*adj spotk*adj Ks*adj];
pchg(t,:)=[dxdt*adj dydt*adj];
loc(t,:)=[Yint Xint];
wchg(t,:)=[dldt*adj dkdt*adj];
som(t,:)=[(Uint*spotk+Vint*spotl+(qxint*spotl-qyint*spotk)/Ks^2)];
stVel(t,:)=[Uint Vint qxint qyint spotk spotl Ks];

if t==5
  run_record=[run_record;kk fr RR];
end
if rem(t,5)==0
alL=[trl nums pchg wchg loc som stVel];alL=double(alL);
%save alL alL -ascii -tabs
if SM==2
eval(sprintf('save jfm200realtrace_wave%d_f%d_root%d alL -ascii -tabs',kk,fr,RR));
else
eval(sprintf('save jfm200realtrace_wave%d_f%d_root%d_nosmooth alL -ascii -tabs',kk,fr,RR));
end
end

end
end
end
end
end
end
end
save raytrace_runrecord5 run_record