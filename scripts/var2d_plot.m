xplot=0:0.75:359.25;
yplot=90:-.75:-90;
[XPLOT,YPLOT]=meshgrid(xplot,yplot(3:239));

figure;
fig1=contour(XPLOT,YPLOT,BetaM(3:239,:));


figure;
fig2=contour(XPLOT,YPLOT,dqbardy(3:239,:));