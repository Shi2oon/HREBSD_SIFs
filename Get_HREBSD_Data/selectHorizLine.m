function [X,Y] = selectHorizLine(x,y,image,V)

Xvec = x(1,:);
Yvec = y(:,1);

H1 = figure; H1=subplot(1,1,1);
hold on
if V.type == 'S'
    imagesc(Xvec,Yvec,image)
%     set(gca,'CLim',[-1.5    1.5]); 
elseif V.type == 'E'
    imagesc(Xvec,Yvec,image)
    set(gca,'CLim',[-5e-3    5e-3]); 
else
    imagesc(Xvec,Yvec,real(log10(image)))
    colormap(jet(256));                 %set(gca,'CLim',[13 15.5]); 
end
axis image;         axis xy;        H1.YDir='reverse';   
colormap jet;       H1.XDir='reverse'; 
% set(gca,'YDir','normal')
axis equal
xlim([min(Xvec) max(Xvec)])
ylim([min(Yvec) max(Yvec)])
colorbar
% set('Name',['Select Horizontal Line (',V.type,'_',num2str(V.i),'_',...
%     num2str(V.j),')'],'NumberTitle','off');
title('Click the Crack start and tip')
xlabel(['x-position [',V.unit,']']);
ylabel(['y-position [',V.unit,']']);
colormap jet

pos = get(gcf,'position');          set(gcf,'position',[100 100 pos(3:4)*2]) 
% X = [24.4412; 35.0460]; Y = [17.3383; 22.4061];
% X = [14.5274; 16.2020 ]; Y = [12.5955; 14.6584];
[X,Y] = ginput(2);
plot(X,Y)
hold off
close all

