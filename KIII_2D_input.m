clc;clear;close all;
addpath(genpath([pwd 'A:\OneDrive - Nexus365\GitHub\3DM\functions']));
set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14)
pname = '';
addpath([pwd '\functions\Get Data'])
% Maps = Calibration_2DKIII(3,1,5); % for calibration
% Maps.results = Maps.Saving;     Dir = Maps;     Dir.Maps = Maps;
[Maps,~] = GetGrainData(pname,'2D_KIII');
[J,KI,KII,KIII] = KIII_2D(Maps);
[~,J,KI,KII,KIII] = Abaqus_2D_KIII(Dir,M4);

%%
clc;clear; addpath(genpath([pwd 'A:\OneDrive - Nexus365\GitHub\3DM\functions']));
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14)
inO = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands';
f{1} = [inO '\S12_1250_R_10nA_20kV_XEBSD_3D\S12_1250_R_10nA_20kV_XEBSD'];
f{4} = [inO '\S2\S2_0750_10nA_20kV_XEBSD'];%minus
f{5} = ['A:\OneDrive - Nexus365\Work\EBSD Data\DIF-GBs\Slips\19-11-06\19-11-06_XEBSD'];%minus
f{6} = [inO '\19-11-07_results\191107_DSS_XEBSD'];
f{7} = 'A:\OneDrive - Nexus365\GitHub\3DM\Calibration\1 WS\1MPa_[um]_XEBSD';

% f{7} = [inO '\S4\S4_1000_10nA_20kV_XEBSD'];
f{8} = [inO '\S4\S4_1200_10nA_20kV_XEBSD'];
f{9} = [inO '\S4\S4_1200_R_10nA_20kV_XEBSD'];

f{10} = [inO '\S5\S5_802um_10nA_20kV_XEBSD'];
f{11} = [inO '\S5\S5_1002um_20nA_10kV_XEBSD'];
f{12} = [inO '\S5\S5_1202um_28nA_20kV_XEBSD'];

f{13} = [inO 'S7\S7_1100_10nA_20kV_XEBSD']; O{13} = '1';
f{14} = [inO 'S7\S7_1100_10nA_20kV_XEBSD']; O{14} = '2';
f{15} = [inO 'S7\S7_1100_10nA_20kV_XEBSD']; O{15} = '3';

f{16} = [inO 'S6\S6_1100_10nA_20kV_XEBSD'];
f{17} = [inO 'S9_1404um_10nA_20kV\200313_InSituDSS_3P_S9_1404um_10nA_20kV_XEBSD'];
f{18} = [inO 'S13_1250_R_10nA_20kV_XEBSD_1_3D\S13_1250_R_10nA_20kV_XEBSD']; O{18} = '1';
f{19} = [inO 'S13_1250_R_10nA_20kV_XEBSD_1_3D\S13_1250_R_10nA_20kV_XEBSD']; O{19} = '2';
f{20} = [inO 'S5_1000_10nA_20kV_XEBSD'];

for iV=7:20
    [Maps,~] = GetGrainData(f{iV},['2D_KIII_' O{iV}]);
    KIII_2D(Maps,Maps.units.S,Maps.units.xy);
end

%%
clc;clear;set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14)
% file  = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S6_1100_10nA_20kV_XEBSD';
% thick = [1,2,3,4,5,7,10,13,17,21,25,29,34,39,44,49,55,66,77,88,100,132,174,233,300,357,450,500];
% file  = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S4_1000um_10nA_20kV_XEBSD';
% thick = [1,3,5,6,7,8,9,13,17,21,25,29,33,44,55,66,77,88,100,122,155,200,250,300,450];
% file = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S9_1404um_10nA_20kV';
% thick = [4,7,11,14,18,22,27,33,39,47,50,56,67,79,94,112,133,159,176,190,232,299,450];%
file = 'A:\OneDrive - Nexus365\Work\EBSD Data\20-11-05 Si Indent\Si_10_2_XEBSD\#WTF';
thick = [17 30 44 56 70 89 101 116 144 167 188 233 256 301 349 406 483 548 670 847 1597];%
zoom=50;
for iV=1:length(thick)
%     if ~exist([file '\Z = ' num2str(thick(iV)) 'nm\Abaqus_2D_KIII.mat'],'file')
        clearvars -except iV thick file Th zoom;
        load([file '\Data2Abaqus_0_reduced'],'Dir')
        Dir.Operation = 'xED';
        Dir.results = [file '\Z = ' num2str(thick(iV)) 'nm'];
        load([Dir.results '\3D_Integrated_Uxy.mat'],'M4');
        [~,J,KI,KII,KIII] = Abaqus_2D_KIII(Dir,M4);
%     end
    load([file '\Z = ' num2str(thick(iV)) 'nm\Abaqus_2D_KIII'])
    %%
    contrs   = 1;
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(2:4))))))+2;
if dic<1;       dic = 1;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:4))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:4)),1),dic);
J.K.true   = round(mean(rmoutliers(J.K.Raw(contrs:4))),dic);
J.K.div    = round(std(rmoutliers(J.K.Raw(contrs:4)),1),dic);

KI.true  = round(mean(rmoutliers(KI.Raw(contrs:4))),dic);
KI.div   = round(std(rmoutliers(KI.Raw(contrs:4)),1),dic);
KII.true = round(mean(rmoutliers(KII.Raw(contrs:4))),dic);
KII.div  = round(std(rmoutliers(KII.Raw(contrs:4)),1),dic);
KIII.true = round(mean(rmoutliers(KIII.Raw(contrs:4))),dic);
KIII.div  = round(std(rmoutliers(KIII.Raw(contrs:4)),1),dic);
    %%
    Th.KI(iV) = KI.true;        Th.KII(iV) = KII.true;
    Th.KIII(iV) = KIII.true;        Th.J(iV) = J.K.true;
    Th.KIv(iV) = KI.div;        Th.KIIv(iV) = KII.div;
    Th.KIIIv(iV) = KIII.div;        Th.Jv(iV) = J.K.div;
end
save([file '\J_KI_II_III_abaqus.mat'],'Th','thick','file');

close all;Kd = [Th.KI(:); Th.KII(:); Th.KIII(:)]; %thick(thick<17)=NaN;
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14)
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;    hold on;
errorbar(thick,Th.KI,Th.KIv,'k--o','MarkerEdgeColor','k','LineWidth',2,'markersize',12);
errorbar(thick,Th.KII,Th.KIIv,'k--s','MarkerEdgeColor','k','LineWidth',...
    2,'MarkerFaceColor','k','markersize',12);
errorbar(thick,Th.KIII,Th.KIIv,'k--<','MarkerEdgeColor','k','LineWidth',2',...
    'markersize',12);
ylabel('K (MPa m^{0.5})'); hold off
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
yyaxis right;set(0,'defaultAxesFontSize',20);
errorbar(thick,Th.J,Th.KIIv,'r--<','MarkerEdgeColor','r','LineWidth',...
    2,'MarkerFaceColor','r','markersize',12);hold on
ylabel('J [J/m^2]');        ylim([0 max(Th.J)+min(Th.J)/4]);
line([zoom zoom],[0 max(Th.J)+min(Th.J)/4],'Color','b','LineStyle','-.','linewidth',3);
hold off
yyaxis left;    hold on;
line([1 1800],[Th.KIII(find(thick==zoom)) Th.KIII(find(thick==zoom))],...
    'Color','b','LineStyle','-.','linewidth',3);
hold off;
xlabel(['Layer thickness [nm]']);        set(gca, 'XScale', 'log'); xlim([0 600])
[~,objh]=legend('K_{I}','K_{II}','K_{III}','@176nm','J_{integral}');%,'location','northoutside','orientation','horizontal');
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 12);
set(gcf,'position',[737 287 955 709]);grid on;  box off;
saveas(gcf, [file '\J_KI_II_III_abaqus.fig']);
saveas(gcf, [file '\J_KI_II_III_abaqus.tif']);close

%%
clc;clear
f{1} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\S6_1250_R_10nA_20kV_XEBSD';
f{2} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\S6_1250_10nA_20kV_XEBSD';
f{3} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\S13_1250_R_10nA_20kV_XEBSD_1';
f{4} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\S13_1250_R_10nA_20kV_XEBSD_2';
f{5} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\S12_1250_R_10nA_20kV_XEBSD_3D';
f{6} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\191106_XEBSD';
f{7} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\x191107_XEBSD';
f{8} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\x1000_10nA_20kV_Slip2';
f{9} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\ExSitu\x1000_10nA_20kV_Slip1';
for iV=7:9
        if ~exist([f{iV} '\Abaqus_2D_KIII.mat'],'file')
        clearvars -except iV f;
        [~,A] = fileparts(f{iV});
        load([f{iV} '\Data2Abaqus_0_reduced'],'Dir')
        Dir.Operation = 'lazy';
        Dir.results = f{iV};
        load([Dir.results '\3D_Integrated_Uxy.mat'],'M4');
        [~,J,KI,KII,KIII] = Abaqus_2D_KIII(Dir,M4);
        end
end
