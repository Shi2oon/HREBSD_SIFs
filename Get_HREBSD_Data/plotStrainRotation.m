function plotStrainRotation(Maps,file)
close
% s1=subplot(3,3,1);  	contourf(Maps.X,Maps.Y,Maps.E11_2,'LineStyle','none');title('\epsilon_{11}');
% axis image; axis off; box off; 
% c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
% s2=subplot(3,3,2);  	contourf(Maps.X,Maps.Y,Maps.E12_2,'LineStyle','none');title('\epsilon_{12}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
% s3=subplot(3,3,3);  	contourf(Maps.X,Maps.Y,Maps.E13_2,'LineStyle','none');title('\epsilon_{13}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
% s5=subplot(3,3,5);  	contourf(Maps.X,Maps.Y,Maps.E22_2,'LineStyle','none');title('\epsilon_{22}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
% s6=subplot(3,3,6);  	contourf(Maps.X,Maps.Y,Maps.E23_2,'LineStyle','none');title('\epsilon_{23}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 
% s9=subplot(3,3,9);  	contourf(Maps.X,Maps.Y,Maps.E33_2,'LineStyle','none');title('\epsilon_{33}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
% addScale([3 3 9],[Maps.X(:) Maps.Y(:)]);
%  
% s4=subplot(3,3,4);  	contourf(Maps.X,Maps.Y,Maps.W12_F1,'LineStyle','none'); 	title('W_{12}');
% axis image; axis off;  box off;
% c  =colorbar;	cU(7,:) = c.Limits;     colorbar off; 
% s7=subplot(3,3,7);  	contourf(Maps.X,Maps.Y,Maps.W13_F1,'LineStyle','none'); 	title('W_{13}');
% axis image; axis off; box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(8,:) = c.Limits;     colorbar off;
% s8=subplot(3,3,8);  	contourf(Maps.X,Maps.Y,Maps.W23_F1,'LineStyle','none'); 	title('W_{23}');
% axis image; axis off;  box off; %set(gca,'Ydir','reverse')
% c  =colorbar;	cU(9,:) = c.Limits;    colorbar off; 
% %
% cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
% caxis(cbax,[min(cU(:)) max(cU(:))]);
% caxis(cbax,[-5e-3 5e-3]);
% h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
% % h.Label.String = '\epsilon'; 
% set([s1 s2 s3 s4 s5 s6 s7 s8 s9],"clim",caxis);% colormap jet;
% %}
% set(gcf,'position',[1 41 1900 900]); 
% saveas(gcf,[file '_Ws and Es.tif'],'tiffn');    
% saveas(gcf,[file '_Ws and Es.fig']);

%%
% X = Maps.X*0.25; Y=Maps.Y*0.25;
% close; Maps.X = unique(Maps.X)*0.25;     Maps.Y = unique(Maps.Y)*0.25;
s1=subplot(3,3,1);  	pcolor(Maps.X,Maps.Y,Maps.E11_2);title('\epsilon_{11}');
axis image; axis off; box off; ;shading interp;
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	pcolor(Maps.X,Maps.Y,Maps.E12_2);title('\epsilon_{12}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;;shading interp;
s3=subplot(3,3,3);  	pcolor(Maps.X,Maps.Y,Maps.E13_2);title('\epsilon_{13}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;;shading interp;
s5=subplot(3,3,5);  	pcolor(Maps.X,Maps.Y,Maps.E22_2);title('\epsilon_{22}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;;shading interp;
s6=subplot(3,3,6);  	pcolor(Maps.X,Maps.Y,Maps.E23_2);title('\epsilon_{23}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; ;shading interp;
s9=subplot(3,3,9);  	pcolor(Maps.X,Maps.Y,Maps.E33_2);title('\epsilon_{33}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;;shading interp;
addScale([3 3 9],[X(:) Y(:)]);
 
s4=subplot(3,3,4);  	pcolor(Maps.X,Maps.Y,Maps.W12_F1); 	title('W_{12}');
axis image; axis off;  box off;
c  =colorbar;	cU(7,:) = c.Limits;     colorbar off; ;shading interp;
s7=subplot(3,3,7);  	pcolor(Maps.X,Maps.Y,Maps.W13_F1); 	title('W_{13}');
axis image; axis off; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(8,:) = c.Limits;     colorbar off;;shading interp;
s8=subplot(3,3,8);  	pcolor(Maps.X,Maps.Y,Maps.W23_F1); 	title('W_{23}');
axis image; axis off;  box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(9,:) = c.Limits;    colorbar off; ;shading interp;
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
% caxis(cbax,[min(cU(:)) max(cU(:))]);
caxis(cbax,[-5e-3 5e-3]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
% h.Label.String = '\epsilon'; 
set([s1 s2 s3 s4 s5 s6 s7 s8 s9],"clim",caxis);% colormap jet;
%}
set(gcf,'position',[1 41 1900 900]); 
saveas(gcf,[file '_W s and Es.tif'],'tiffn');    
saveas(gcf,[file '_W s and Es.fig']);    