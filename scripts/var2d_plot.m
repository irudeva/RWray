xplot=0:0.75:359.25;
yplot=90:-.75:-90;
[XPLOT,YPLOT]=meshgrid(xplot,yplot(3:239));

fig1=contour(XPLOT,YPLOT,BetaM(3:239,:));
title('BetaM')

figure;
fig2=contour(XPLOT,YPLOT,dqbardy(3:239,:));
title('qbardy')

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