xplot=0:0.75:359.25;
yplot=90:-.75:-90;
[XPLOT,YPLOT]=meshgrid(xplot,yplot(3:239));
[XPLOT1,YPLOT1]=meshgrid(xplot,yplot(4:238));

%fig1=contour(XPLOT,YPLOT,BetaM(3:239,:));
%title('BetaM')

%figure;
%fig2=contour(XPLOT1,YPLOT1,dqbardx(4:238,:));
%title('qbardx')

%figure;
%fig3=contour(XPLOT1,YPLOT1,qbar(4:238,:));
%title('qbar')

figure;
fig=contour(XPLOT1,YPLOT1,u0(4:238,:));
title('u JJA')

%v=[-1:0.5:1];
%v=v*1.e-14;
%fig3=contour(XPLOT,YPLOT,qbar(3:239,:),v);
%title('qbardy')


%figure;
%v=[-1.5:0.25:1.5];
%v=v*1.e-4;
%fig3=contour(XPLOT,YPLOT,qbar(3:239,:),v);
%title('qbar');

%figure;
%v=[-1.5:0.25:1.5];
%v=v*1.e-4;
%fig4=contour(XPLOT,YPLOT,d2qbardxy(3:239,:),v);
%title('d2qbardxy');