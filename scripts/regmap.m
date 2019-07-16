addpath /Applications/MATLAB_R2016a.app/m_map/

% reg = ['C_Europe ';'E_Europe ';'W_Siberia';'E_Siberia';'W_China  ';'E_China  ';'W_America';'E_America' ];
reg = ['C Europe ';'E Europe ';'W Siberia';'E Siberia';'W China  ';'E China  ';'W America';'E America';'Kara Bar ';'ESib Chuk' ];

lonreg1 = [ 5, 30, 70, 100,  80, 100, 230, 260, 30, 160 ];
lonreg2 = [25, 50, 90, 120, 100, 120, 255, 280, 70, 200 ];
latreg1 = [40, 45, 50,  55,  25,  22,  45,  35, 70,  70 ];
latreg2 = [55, 60, 65,  70,  45,  42,  60,  50, 80,  80 ];


clf;
axes('position',[.1 .2 .8 .7]);
m_proj('mercator','lon',[0 360],'lat',[0 85]);

for ireg = 1:1:length(reg(:,1))

bndry_lon=[lonreg1(ireg) lonreg1(ireg) lonreg2(ireg) lonreg2(ireg) lonreg1(ireg)];
bndry_lat=[latreg1(ireg) latreg2(ireg) latreg2(ireg) latreg1(ireg) latreg1(ireg)];



%         % [cs,h]= m_contourf(lon,lat(6:176),(real(sqrt(rdivide(sBetaM1,sUbarM1)))*rad));
%         [cs,h]= m_contourf(lon,lat,Ks,[0:1:20]);

%         hold on;

         m_coast('linewidth',1,'color','r');
         m_grid('linewi',2,'linest','none','tickdir','out','fontsize',12);;

         m_line(bndry_lon,bndry_lat,'linewi',2,'color','k');     % Area outline ...
%          m_text((lonreg1(ireg)+lonreg2(ireg))/2 ,(latreg1(ireg)+latreg2(ireg))/2,5,{'Pacific','Ocean'},'fontsize',18);
         m_text(lonreg1(ireg)+1 ,(latreg1(ireg)+latreg2(ireg))/2,5,reg(ireg,:),'fontsize',12);

%          title(sprintf('Ks, %s%d',bgf,fyr1),'fontsize',14);
%          ax=m_contfbar([.2 .8],-.1,cs,h);
%          set(ax,'fontsize',12)

end

% %         figout     = sprintf('../output/Ks/Ks.%s%d-%d.png',bgf,fyr1,lyr); 
%         figout     = sprintf('../output/Ks/Ks.%s%d.png',bgf,fyr1); 
