function plotStressandRot(S)
S.S11(isnan(S.GND))= NaN;   S.S22(isnan(S.GND))= NaN;   S.W21(isnan(S.GND))= NaN;   
S.W31(isnan(S.GND))= NaN;   S.W32(isnan(S.GND))= NaN;   S.S13(isnan(S.GND))= NaN;
S.S23(isnan(S.GND))= NaN;   close all;                  S.S12(isnan(S.GND))= NaN; 
O.X = S.X; O.Y = S.Y;
close; S.X = unique(S.X);     S.Y = unique(S.Y);

%%
set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14)
s1=subplot(3,3,1);  	pcolor(S.X,S.Y,S.S11);
title('\sigma_{xx}','fontsize',19);set(gca,'Ydir','reverse')
axis image; axis off; colormap jet;     box off; shading interp;
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	pcolor(S.X,S.Y,S.S12);
title('\sigma_{xy}','fontsize',19); shading interp;
axis image; axis off; colormap jet;     box off; 
set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	pcolor(S.X,S.Y,S.S13);
title('\sigma_{xz}','fontsize',19); shading interp;
axis image; axis off; colormap jet;     box off; 
set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	pcolor(S.X,S.Y,S.S22);
title('\sigma_{yy}','fontsize',19); shading interp;
axis image; axis off; colormap jet;     box off; 
set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	pcolor(S.X,S.Y,S.S23);
title('\sigma_{yz}','fontsize',19); shading interp;
axis image; axis off; colormap jet;     box off; 
set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 

s9=subplot(3,3,9);  	pcolor(S.X,S.Y,log10(S.GND)); 	
title('\rho_{GNDs}','fontsize',19); axis image; axis off; colormap jet; box off; 
set(gca,'Ydir','reverse'); shading interp;
addScale([3 3 9],[O.X(:) O.Y(:)]);
subplot(3,3,9); c=colorbar;  c.Label.String = ['\rho_{G} [log_{10} m^{-2}]'];
c.Position = [0.8774 0.1135 0.0111 0.1942];
 
s4=subplot(3,3,4);  	pcolor(S.X,S.Y,S.W21); 	
title('\omega_{yx}','fontsize',19);set(gca,'Ydir','reverse')
axis image; axis off;  box off; colormap jet; shading interp;
c  =colorbar;	cu(1,:) = c.Limits;     colorbar off; 
s7=subplot(3,3,7);  	pcolor(S.X,S.Y,S.W31); 	 shading interp;
title('\omega_{zx}','fontsize',19);set(gca,'Ydir','reverse')
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cu(2,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	   
pcolor(S.X,S.Y,S.W32);title('\omega_{zy}','fontsize',19);
axis image; axis off; colormap jet; box off;  shading interp;
set(gca,'Ydir','reverse')
c  =colorbar;	cu(3,:) = c.Limits;     colorbar off;

cu(abs(cu)==1)=0; set([s4 s7 s8],"clim",[-0.05 0.05]); 
subplot(3,3,4); c=colorbar;  c.Label.String = [ '\omega [Rad]']; 
c.Position = [00.1236 0.1133 0.0112 0.4611];

cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[-1.5 1.5]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.8774 0.4178 0.0111 0.4800] );
h.Label.String = '\sigma [GPa]'; set([s1 s2 s3 s5 s6],"clim",caxis); 
set(gcf,'position',[1 41 1900 900]); 