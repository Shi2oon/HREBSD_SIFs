clc;clear
Map = Calibration_2DKIII(3,1,5);
count=0;
for iX=0.001:0.0001:0.1
    Maps = Map;
    for iV=1:3
        for iO=1:3
            eval(sprintf('a=Map.E%d%d;',iV,iO));% Define our signal - a column vector.
            randomNoise = randn(size(a)) .* iX * a;% Create the noise values that we'll add to a.
            % Add noise to a to make an output column vector.
            eval(sprintf('Maps.E%d%d = a + randomNoise;',iV,iO));
        end
    end
    subplot(1,2,1); imagesc(Map.E11); axis image
    subplot(1,2,2); imagesc(Maps.E11); axis image
    alldata = [Maps.X(:) Maps.Y(:) Maps.Z(:) Maps.E11(:) Maps.E12(:) Maps.E13(:)...
        Maps.E21(:) Maps.E22(:) Maps.E23(:) Maps.E31(:) Maps.E32(:) Maps.E33(:)];
    [J,KI,KII,KIII] = KIII_2D(alldata,Maps); % as desigignated maps
    count = count+1;
    j(count) = J.true;           je(count) = J.div;
    ki(count) = KI.true;          kie(count) = KI.div;
    kii(count) = KII.true;         kiie(count) = KII.div;
    kiii(count) = KIII.true;        kiiie(count) = KIII.div;
end

%%
Ki = [0.001 1.2 3.1 3.9 5.6 6.2 7.1 8.1 8.7 9 9.4];
Kii = [0.001 1.8 4.8 5.9 6.8 8.1 9.4];
Kiii = [1 4 5.2 6.6 7.9 9.3];

for iV=1:length(Ki)
    [~,ind_Ki(iV)] = min(abs(Contour-Ki(iV)));
end
%%
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14)
Contour = (0.001:0.0001:0.1)*100; close all

ft = fittype( 'poly3' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf 0];
opts.Robust = 'LAR';
opts.Upper = [Inf Inf 0];


f1=fit([Contour'],[abs(ki'./ki(1)-1)*100],ft);
f1e=fit([Contour'],100-[abs(kie'./ki(1)-1)*100],ft);
p1=plot(f1);hold on;
p1.LineWidth = 4;       p1.Color = 'k';
f2=fit([Contour'],[abs(kii'./kii(1)-1)*100],ft);
f2e=fit([Contour'],100-[abs(kiie'./kii(1)-1)*100],ft);
p1=plot(f2);hold on;
p1.LineWidth = 4;       p1.Color = 'b';
f3=fit([Contour'],[abs(kiii'./kiii(1)-1)*100],ft);
f3e=fit([Contour'],100-[abs(kiiie'./kiii(1)-1)*100],ft);
p1=plot(f3);hold on;
p1.LineWidth = 4;       p1.Color = 'c';
f4=fit([Contour'],[abs(j'./j(1)-1)*100],ft);
f4e=fit([Contour'],100-[abs(je'./j(1)-1)*100],ft);
p1=plot(f4);
p1.LineWidth = 4;       p1.Color = 'r';

p1=plot(f1(1:10)+f1e(1:10));
p1.LineWidth = 2;       p1.Color = 'k'; p1.LineStyle = '--';
p1=plot(f2(1:10)+f2e(1:10));
p1.LineWidth = 2;       p1.Color = 'b'; p1.LineStyle = '--';
p1=plot(f3(1:10)+f3e(1:10));
p1.LineWidth = 2;       p1.Color = 'c'; p1.LineStyle = '--';
p1=plot(f4(1:10)+f4e(1:10));
p1.LineWidth = 2;       p1.Color = 'r'; p1.LineStyle = '--';hold off;
% plot(Contour,abs(kii./kii(1)-1)*100,'m.');
% plot(Contour,abs(kiii./kiii(1)-1)*100,'b.');
% plot(Contour,abs(j./j(1)-1)*100,'k.');
ylabel('%Error');        axis tight; %ylim([0 30])

% stabi = mean([abs(je./j(1)-1)*100;      abs(kie./ki(1)-1)*100;  ...
%       abs(kiie./kii(1)-1)*100;  abs(kiiie./kiii(1)-1)*100]);
% plot(Contour,abs(ki./ki(1)-1),kie./ki(1),'r.');hold on;

grid off; xlim([0   10]); ylim([0 60]); yticks([0:10:90]); xticks([0:1:10])
[~,objh]=legend('K_{I}','K_{II}','K_{III}','J_{integral}','Upper Bound','location','northwest');
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 28);
xlabel('Random noise %'); box off
set(gcf,'position',[203 196 556 642]);box off;
saveas(gcf,'A:\OneDrive - Nexus365\GitHub\3DM\functions\2D_KIII\Crack_tip_4.tif'); 
saveas(gcf,'A:\OneDrive - Nexus365\GitHub\3DM\functions\2D_KIII\Crack_tip_4.fig'); close 

