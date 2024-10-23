function [K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(alldata,MatProp)
close all;

% This code decompose the Stress intesity factors from strain maps
% directly without the need for integration
% Start with generating strain data using the calibration code and
% then use the output for the main function
% The code is self contained and does not need extra functions

% for this function to work properly, spacing between potins in x and y
% should be the same and the crack should be at the centre (this can be
% done inside this fucntion also)


% this functions accept data from HR-EBSD as formuated with 'GetData'
% function. The crack needs to be on the centre exactly for the code to
% work and give good results. This code already include the assumption of
% sigma 33 free == 0
% stress in Pa, E in Pa, distance in m
% use variable 'Stiffness' if you are using anistropic material or E and nu
% for istroropic material
% need some unit clbartion specially for mm

% Another option is to input the defromation gradient components as a vector matrix with
% 9 columns the first three columns are the x, y and z coordinate in meters. z
% cooridnate can be a zero column. the 4th to the 9th column are the strain
% components arranged as
% Maps = [X(:) Y(:) Z(:) E11(:) E12(:) E13(:) E21(:) E22(:) E23(:) E31(:) E32(:) E33(:)];
% or A for deformation gradient tensor

% if the map is a 2D defromation gradient map then zero all out of the plane components

% the material paramters as MatProp.E for Young's Modulus and MatProp.nu
% for Possions ratio. or as a stifness matrix all in Pa
% xEBSD will assume and solve for  plane strain conditions with sigma33 == 0 + no
% volumetric change.

% Example:
% [MatProp,~,alldata] = Calibration_2DKIII(3,1,2);
% [J,KI,KII,KIII] = KIII_2D(alldata,MatProp);% or just MatProp
% [J,KI,KII,KIII] = KIII_2D(MatProp); % as desigignated maps

% look to the strrcutrure of each variable to see the differance

if isfield(MatProp,'Operation')
    if strcmpi(MatProp.Operation, 'DIC')
        [~,RawData ] = reshapeData(alldata);
        stepsize = unique(round(diff(unique(RawData.Y1(:))),4));
        [RawData.E11,RawData.E12,RawData.E13] = crackgradient(RawData.Ux,stepsize);
        [RawData.E21,RawData.E22,RawData.E23] = crackgradient(RawData.Uy,stepsize);
        U (:,:,1) = RawData.Ux;        U (:,:,2) = RawData.Uy;

        if size(alldata,2)==6
            [RawData.E31,RawData.E32,RawData.E33] = crackgradient(RawData.Uz,stepsize);

            alldata = [RawData.X1(:) RawData.Y1(:) zeros(size(RawData.Y1(:))) ...
                RawData.E11(:) RawData.E12(:) RawData.E13(:) ...
                RawData.E21(:) RawData.E22(:) RawData.E23(:) ...
                RawData.E31(:) RawData.E32(:) RawData.E33(:)];
            U (:,:,3) = RawData.Uz;
        else
            alldata = [RawData.X1(:) RawData.Y1(:) zeros(size(RawData.Y1(:))) ...
                RawData.E11(:) RawData.E12(:) zeros(size(RawData.Y1(:)))...
                RawData.E21(:) RawData.E22(:) zeros(size(RawData.Y1(:))) ...
                zeros(size(RawData.Y1(:))) zeros(size(RawData.Y1(:))) zeros(size(RawData.Y1(:)))];
        end

    end
end

if size(alldata,2) == 5 % 2D deformation gradient data
    alldata = [alldata(:,1) alldata(:,2) zeros(size(alldata(:,2))) ...
        alldata(:,3) alldata(:,5) zeros(size(alldata(:,2))) ...
        alldata(:,5) alldata(:,4) zeros(size(alldata(:,2))) ...
        zeros(size(alldata(:,2))) zeros(size(alldata(:,2))) ...
        zeros(size(alldata(:,2)))];
end
[~,Maps]=reshapeDefromationGradient(alldata);
if size(MatProp,1) == 6
    Maps.Stiffness = MatProp;
else
    if isfield(MatProp,'Stiffness')
        Maps.Stiffness = MatProp.Stiffness;
    else
        Maps.E  = MatProp.E;
        Maps.nu = MatProp.nu;
    end
    Maps.stressstat = MatProp.stressstat;
    Maps.units.xy = MatProp.units.xy;
    Maps.units.St = MatProp.units.St;
    if isfield(MatProp,'SavingD')
        Maps.SavingD = MatProp.SavingD;
    end
end
%
%% prepare Data
imagesc(Maps.du11);
axis tight; axis image; axis off; colormap jet
set(gcf,'position',[737 287 955 709]);
%
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'N';         % Use the TeX interpreter to format the question
quest            = 'Do you want to Crop and Centre the Crack tip';
answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if strcmpi(answer,'Y') % crop data
    [Crop] = CroppingEqually(Maps);
    Maps.X    = Crop.X;      Maps.Y   = Crop.Y;        Maps.Z   = Crop.Z;
    Maps.du11 = Crop.du11;  Maps.du12 = Crop.du12;    Maps.du13 = Crop.du13;
    Maps.du21 = Crop.du21;  Maps.du22 = Crop.du22;    Maps.du23 = Crop.du23;
    Maps.du31 = Crop.du31;  Maps.du32 = Crop.du32;    Maps.du33 = Crop.du33;
end
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'L';         % Use the TeX interpreter to format the question
quest            = 'Is the crack on your left or right ?';
answer           = questdlg(quest,'Boundary Condition','L','R', opts);
if strcmpi(answer,'R') % crop data
    %}
    Maps.du11 = flip(flip(Maps.du11,1),2);    Maps.du12 = flip(flip(Maps.du12,1),2);
    Maps.du13 = flip(flip(Maps.du13,1),2);
    Maps.du21 = flip(flip(Maps.du21,1),2);    Maps.du22 = flip(flip(Maps.du22,1),2);
    Maps.du23 = flip(flip(Maps.du23,1),2);
    Maps.du31 = flip(flip(Maps.du31,1),2);    Maps.du32 = flip(flip(Maps.du32,1),2);
    Maps.du33 = flip(flip(Maps.du33,1),2);
end
close
%}
%%
switch Maps.units.St
    case 'Pa'
        Saf = 1;
    case 'KPa'
        Saf = 1e3;
    case 'MPa'
        Saf = 1e6;
    case 'GPa'
        Saf = 1e9;
end

if isfield(Maps,'Stiffness')
    Maps.Stiffness = Maps.Stiffness*Saf;
    [Maps.E,Maps.nu,Maps.G,Maps.Co] = effectiveE_nu(Maps.Stiffness); % in Pa
else
    Maps.E = Maps.E*Saf;
    Maps.G = Maps.E/(2*(1 + Maps.nu));
    if strcmpi(Maps.stressstat,'plane_strain')
        Maps.E = Maps.E/(1-Maps.nu^2);% for HR-EBSD plane stress conditions
    end
end
switch Maps.units.xy
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e-3;
    case 'um'
        saf = 1e-6;
    case 'nm'
        saf = 1e-9;
end
Maps.X = Maps.X*saf;
Maps.Y = Maps.Y*saf;

Maps.stepsize = unique(round(diff(unique(Maps.Y(:))),4));
Maps.units.St = 'Pa';        Maps.units.xy = 'um';
DataSize = [size(Maps.du11),1];

%% Decomposition method.
[du_dx,E,S] = decomposeDU(Maps);
Wd = 0.5*(E(:,:,1,1,:).*S(:,:,1,1,:) + E(:,:,1,2,:).*S(:,:,1,2,:) + E(:,:,1,3,:).*S(:,:,1,3,:)...
    + E(:,:,2,1,:).*S(:,:,2,1,:) + E(:,:,2,2,:).*S(:,:,2,2,:) + E(:,:,2,3,:).*S(:,:,2,3,:)...
    + E(:,:,3,1,:).*S(:,:,3,1,:) + E(:,:,3,2,:).*S(:,:,3,2,:) + E(:,:,3,3,:).*S(:,:,3,3,:));
Wd = squeeze(Wd);
%{
% Decomposed Plots
plot_DecomposeddU(du_dx,Maps);
if isfield(Maps,'SavingD')
    saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_du.fig']);
    saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_du.tif']);  close
end

if isfield(Maps,'E11')
else
    Maps.E11 = Maps.du11; Maps.E12 = Maps.du12;
    Maps.E31 = Maps.du31; Maps.E22 = Maps.du22;
    Maps.E32 = Maps.du32; Maps.E33 = Maps.du33;
end
plot_DecomposedStrain(E(:,:,1,1,:),E(:,:,2,2,:),E(:,:,3,3,:),E(:,:,1,2,:),...
    E(:,:,1,3,:),E(:,:,2,3,:),Maps);
if isfield(Maps,'SavingD')
    saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_Strain.fig']);
    saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_Strain.tif']);  close
end

if isfield(Maps,'S11')
    plot_DecomposedStess(S(:,:,1,1,:),S(:,:,2,2,:),S(:,:,3,3,:),S(:,:,1,2,:),...
        S(:,:,1,3,:),S(:,:,2,3,:),Maps,Saf);
    if isfield(Maps,'SavingD')
        saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_Stress.fig']);
        saveas(gcf, [fileparts(Maps.SavingD) '\Decomposed_Stress.tif']);  close
    end
end
%}
%%
% An approach to calculate the J -integral by digital image correlation
% displacement field measurement, FFEMS (2012), 35, ?971-984
% Displacement gradient
% Generate q field
celw = 1;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(Maps.stepsize*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*Maps.stepsize^2;
dq1d(1,:,:) = dQdX;  dq1d(2,:,:) = dQdY;
Am = [1,0; 0, 1];

for jj = 1:2
    for ii=1:3
        termJ(:,:,:,ii) = squeeze(S(:,:,ii,jj,:)).*squeeze(du_dx(:,:,ii,1,:));
    end
    JAdr(:,:,:,jj) =(sum(termJ,4) -Wd.*Am(1,jj)).*squeeze(dq1d(jj,:,:)).*dA;
end
JAd = sum(JAdr,4);
%}
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2));
    J.Raw(:,ii,:) = sum(sum(JAd.*areaID,'omitnan'),'omitnan');
end

%% Equivalent SIF
% to avoid imaginary number (needs to be solved so the code could work for
% compressive fields
J.Raw = J.Raw;
J.KRaw(1:2,:) = sqrt(J.Raw(1:2,:)*Maps.E);
J.KRaw(3,:) = sqrt(J.Raw(3,:)*2*Maps.G);      % Mode III
J.JRaw = J.Raw;
J.Raw = sum(J.Raw);
%
%%
figure; plot(J.Raw); legend('J')%trim acess
set(gcf,'position',[98 311 1481 667])
text(1:length(J.Raw),J.Raw,string([1:length(J.Raw)]))
%}
pause(0.1)
oh = input('where to cut the contour? '); close

%%
J.Raw    = J.Raw(1:oh);
J.KRaw   = J.KRaw(:,1:oh);
KI.Raw   = J.KRaw(1,:)*1e-6;
KII.Raw  = J.KRaw(2,:)*1e-6;
KIII.Raw = J.KRaw(3,:)*1e-6;

%%
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<2;       dic = 2;    end
J.true   = round(mean(J.Raw(contrs:end)),dic);
J.div    = round(std(J.Raw(contrs:end),1),dic);
KI.true  = round(mean(KI.Raw(contrs:end)),dic);
KI.div   = round(std(KI.Raw(contrs:end),1),dic);
KII.true = round(mean(KII.Raw(contrs:end)),dic);
KII.div  = round(std(KII.Raw(contrs:end),1),dic);
KIII.true= round(mean(KIII.Raw(contrs:end)),dic);
KIII.div = round(std(KIII.Raw(contrs:end),1),dic);
K.Raw    = sqrt(J.Raw*Maps.E)*1e-6;
K.true   = round(mean(K.Raw(contrs:end)),dic);
K.div    = round(std(K.Raw(contrs:end),1),dic);
%
plot_JKIII(KI,KII,KIII,J,Maps.stepsize/saf,Maps.units.xy)
if isfield(Maps,'SavingD')
    saveas(gcf, [fileparts(Maps.SavingD) '\J_K.fig']);
    saveas(gcf, [fileparts(Maps.SavingD) '\J_K.tif']);  close all
    save([fileparts(Maps.SavingD) '\KIII_2D.mat'],'J','K','KI','KII','KIII','Maps');
end
%}

%% undecomposed M and J integrals
%% Wd.
clear JAdr MAdr E S du_dx Wd
[du_dx,E,S] = DU(Maps);
Wd = 0.5*(E(:,:,1,1).*S(:,:,1,1) + E(:,:,1,2).*S(:,:,1,2) + E(:,:,1,3).*S(:,:,1,3)...
    + E(:,:,2,1).*S(:,:,2,1) + E(:,:,2,2).*S(:,:,2,2) + E(:,:,2,3).*S(:,:,2,3)...
    + E(:,:,3,1).*S(:,:,3,1) + E(:,:,3,2).*S(:,:,3,2) + E(:,:,3,3).*S(:,:,3,3));
Wd = squeeze(Wd);
%}
%% M-integral
% displacement field measurement, FFEMS (2012), 35, ?971-984
% Displacement gradient
xx = [Maps.stepsize Maps.stepsize];
for kk = 1:2
    for jj = 1:2
        clear termJ termM
        for ii=1:3
            termJ(:,:,ii) = squeeze(S(:,:,ii,jj)).*squeeze(du_dx(:,:,ii,kk));
            termM(:,:,ii) = squeeze(S(:,:,ii,jj)).*squeeze(du_dx(:,:,ii,kk)).*xx(kk);
        end
        JAdr(:,:,kk,jj) =(sum(termJ,3) -Wd.*Am(1,jj)).*squeeze(dq1d(jj,:,:)).*dA;
        MAdr(:,:,kk,jj) = (sum(termM,3)-Wd.*xx(jj).*Am(kk,jj)).*squeeze(dq1d(jj,:,:)).*dA;
    end
end

JAd = sum(JAdr,4);
MAd = sum(MAdr,4); %J/m
%}
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2));
    M.Raw(:,ii,:) = sum(sum(MAd.*areaID,'omitnan'),'omitnan');
    J.vectorial(:,ii,:) = sum(sum(JAd.*areaID,'omitnan'),'omitnan');
end

% to avoid imaginary number (needs to be solved so the code could work for
% compressive fields
J.vectorial = J.vectorial(:,1:oh);
J.vectorial_true   = round(mean(J.vectorial(:,contrs:end),2),dic);
J.vectorial_div    = round(std(J.vectorial(:,contrs:end),1,2),dic);

dic = real(ceil(-log10(nanmean(rmoutliers(M.Raw(contrs:end))))))+2;
if dic<2;       dic = 2;    end
M.Raw = M.Raw(:,1:oh);
M.true  = round(mean(M.Raw(:,contrs:end),2),dic);
M.div   = round(std(M.Raw(:,contrs:end),1,2),dic);
%
plot_JM(M,J,Maps.stepsize/saf,Maps.units.xy,'M')
if isfield(Maps,'SavingD')
    saveas(gcf, [fileparts(Maps.SavingD) '\J_M.fig']);
    saveas(gcf, [fileparts(Maps.SavingD) '\J_M.tif']);  close all
    save([fileparts(Maps.SavingD) '\KIII_2D.mat'],'J','K','KI','KII','KIII','Maps','M');
end
%}

%% L-integral (not correct!, see https://doi.org/10.1007/s00707-014-1152-y
%{
if exist('U','var')
    emij=[0,1;-1,0]; %% for
    for ll=1:2
        for kk = 1:2
            for jj = 1:2
                clear termL
                for ii=1:3
                    termL1(:,:,ii) = squeeze(S(:,:,ii,ll)).*squeeze(du_dx(:,:,ii,kk)).*xx(jj);
                end
                LAdr(:,:,kk,ll,jj)=(Wd.*xx(jj).*Am(kk,ll) + ...
                    squeeze(S(:,:,kk,ll)).*squeeze(U(:,:,jj)).*saf-sum(termL1,3))...
                    .*squeeze(dq1d(ll,:,:)).*dA.*emij(kk,jj);
            end
        end
    end
    LAd = sum(sum(sum(LAdr,6),5),4); %J/m
    for ii = 1:mid-1
        areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2));
        L.Raw(:,ii,:) = sum(sum(LAd.*areaID,'omitnan'),'omitnan');
    end

dic = real(ceil(-log10(nanmean(rmoutliers(L.Raw(contrs:end))))))+2;
if dic<2;       dic = 2;    end
L.Raw = abs(L.Raw(:,1:oh));
L.true  = round(mean(L.Raw(:,contrs:end),2),dic);
L.div   = round(std(L.Raw(:,contrs:end),1,2),dic);

plot_JM(L,J,Maps.stepsize/saf,Maps.units.xy,'L')
end
%}
end


%% sub-functions
function [du_dx,De_E,De_S] = decomposeDU(Maps)
for iV=1:3
    for xi=1:3
        eval(sprintf('A(:,:,iV,xi) = Maps.du%d%d;',iV,xi));
    end
end

%% decompostion
% Mode I
du_dx(:,:,1,1,1) = 0.5*(squeeze(A(:,:,1,1)) + flipud(squeeze(A(:,:,1,1))));
du_dx(:,:,1,2,1) = 0.5*(squeeze(A(:,:,1,2)) - flipud(squeeze(A(:,:,1,2))));
du_dx(:,:,1,3,1) = 0.5*(squeeze(A(:,:,1,3)) + flipud(squeeze(A(:,:,1,3))));

du_dx(:,:,2,1,1) = 0.5*(squeeze(A(:,:,2,1)) - flipud(squeeze(A(:,:,2,1))));
du_dx(:,:,2,2,1) = 0.5*(squeeze(A(:,:,2,2)) + flipud(squeeze(A(:,:,2,2))));
du_dx(:,:,2,3,1) = 0.5*(squeeze(A(:,:,2,3)) - flipud(squeeze(A(:,:,2,3))));

du_dx(:,:,3,1,1) = 0.5*(squeeze(A(:,:,3,1)) + flipud(squeeze(A(:,:,3,1))));
du_dx(:,:,3,2,1) = 0.5*(squeeze(A(:,:,3,2)) - flipud(squeeze(A(:,:,3,2))));
du_dx(:,:,3,3,1) = 0.5*(squeeze(A(:,:,3,3)) + flipud(squeeze(A(:,:,3,3))));

% Mode II
du_dx(:,:,1,1,2) = 0.5*(squeeze(A(:,:,1,1)) - flipud(squeeze(A(:,:,1,1))));
du_dx(:,:,1,2,2) = 0.5*(squeeze(A(:,:,1,2)) + flipud(squeeze(A(:,:,1,2))));
du_dx(:,:,1,3,2) = zeros(size(squeeze(A(:,:,1,1))));

du_dx(:,:,2,1,2) = 0.5*(squeeze(A(:,:,2,1)) + flipud(squeeze(A(:,:,2,1))));
du_dx(:,:,2,2,2) = 0.5*(squeeze(A(:,:,2,2)) - flipud(squeeze(A(:,:,2,2))));
du_dx(:,:,2,3,2) = zeros(size(squeeze(A(:,:,1,1))));

du_dx(:,:,3,1,2) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,3,2,2) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,3,3,2) = 0.5*(squeeze(A(:,:,3,3)) - flipud(squeeze(A(:,:,3,3))));

% Mode III
du_dx(:,:,1,1,3) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,1,2,3) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,1,3,3) = 0.5*(squeeze(A(:,:,1,3)) - flipud(squeeze(A(:,:,1,3))));

du_dx(:,:,2,1,3) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,2,2,3) = zeros(size(squeeze(A(:,:,1,1))));
du_dx(:,:,2,3,3) = 0.5*(squeeze(A(:,:,2,3)) + flipud(squeeze(A(:,:,2,3))));

du_dx(:,:,3,1,3) = 0.5*(squeeze(A(:,:,3,1)) - flipud(squeeze(A(:,:,3,1))));
du_dx(:,:,3,2,3) = 0.5*(squeeze(A(:,:,3,2)) + flipud(squeeze(A(:,:,3,2))));
du_dx(:,:,3,3,3) = zeros(size(squeeze(A(:,:,1,1))));
%
% strain from displacement gradient
for i=1:3
    for j=1:3
        for M=1:3
            De_E(:,:,i,j,M) = 0.5*(du_dx(:,:,i,j,M)+du_dx(:,:,j,i,M));
        end
    end
end

%%
if ~isfield(Maps,'Stiffness') % linear istropic material
    % Chauchy stress tensor, assming linear-elastic, isotropic material for
    De_S(:,:,1,1,:) = Maps.E/(1-Maps.nu^2)*(De_E(:,:,1,1,:) + ...
        Maps.nu*(De_E(:,:,2,2,:) + De_E(:,:,3,3,:)));
    De_S(:,:,2,2,:) = Maps.E/(1-Maps.nu^2)*(De_E(:,:,2,2,:) + ...
        Maps.nu*(De_E(:,:,1,1,:) + De_E(:,:,3,3,:)));
    De_S(:,:,3,3,:) = Maps.E/(1-Maps.nu^2)*(De_E(:,:,3,3,:) + ...
        Maps.nu*(De_E(:,:,1,1,:) + De_E(:,:,2,2,:)));
    De_S(:,:,1,2,:) = Maps.G*(De_E(:,:,1,2,:) + De_E(:,:,2,1,:));
    De_S(:,:,2,1,:) = Maps.G*(De_E(:,:,1,2,:) + De_E(:,:,2,1,:));
    De_S(:,:,1,3,:) = Maps.G*(De_E(:,:,1,3,:) + De_E(:,:,3,1,:));
    De_S(:,:,3,1,:) = Maps.G*(De_E(:,:,1,3,:) + De_E(:,:,3,1,:));
    De_S(:,:,2,3,:) = Maps.G*(De_E(:,:,2,3,:) + De_E(:,:,3,2,:));
    De_S(:,:,3,2,:) = Maps.G*(De_E(:,:,2,3,:) + De_E(:,:,3,2,:));
else% ansitropic material
    for yi = 1:size(A,1)
        for xi = 1:size(A,2)
            for iV = 1:3
                e_voight = [De_E(yi,xi,1,1,iV);  De_E(yi,xi,2,2,iV);...
                    De_E(yi,xi,3,3,iV);   ...
                    De_E(yi,xi,1,2,iV) + De_E(yi,xi,2,1,iV); ...
                    De_E(yi,xi,1,3,iV) + De_E(yi,xi,3,1,iV); ...
                    De_E(yi,xi,2,3,iV) + De_E(yi,xi,3,2,iV)];
                s_voight = Maps.Stiffness*e_voight;
                De_S(yi,xi,:,:,iV) = [s_voight(1),s_voight(4),s_voight(5);
                    s_voight(4),s_voight(2),s_voight(6);
                    s_voight(5),s_voight(6),s_voight(3)];
                De_E(yi,xi,:,:,iV)=[e_voight(1),  e_voight(4)/2,e_voight(5)/2;
                    e_voight(4)/2,e_voight(2),  e_voight(6)/2;
                    e_voight(5)/2,e_voight(6)/2,e_voight(3)];
            end
        end
    end
end
end

%%
function [A, De_E, De_S] = DU(Maps)
for iV = 1:3
    for xi = 1:3
        eval(sprintf('A(:,:,iV,xi) = Maps.du%d%d;', iV, xi));  % Assuming Maps.E11, Maps.E12, etc.
    end
end

% Strain from displacement gradient
for i = 1:3
    for j = 1:3
        De_E(:,:,i,j) = 0.5 * (A(:,:,i,j) + A(:,:,j,i));
    end
end

%% Stress Tensor Calculation
if ~isfield(Maps, 'Stiffness') % linear isotropic material
    % Cauchy stress tensor, assuming linear-elastic, isotropic material
    De_S(:,:,1,1) = Maps.E / (1 - Maps.nu^2) * (De_E(:,:,1,1) + Maps.nu * (De_E(:,:,2,2) + De_E(:,:,3,3)));
    De_S(:,:,2,2) = Maps.E / (1 - Maps.nu^2) * (De_E(:,:,2,2) + Maps.nu * (De_E(:,:,1,1) + De_E(:,:,3,3)));
    De_S(:,:,3,3) = Maps.E / (1 - Maps.nu^2) * (De_E(:,:,3,3) + Maps.nu * (De_E(:,:,1,1) + De_E(:,:,2,2)));
    De_S(:,:,1,2) = Maps.G * (De_E(:,:,1,2) + De_E(:,:,2,1));
    De_S(:,:,2,1) = De_S(:,:,1,2);
    De_S(:,:,1,3) = Maps.G * (De_E(:,:,1,3) + De_E(:,:,3,1));
    De_S(:,:,3,1) = De_S(:,:,1,3);
    De_S(:,:,2,3) = Maps.G * (De_E(:,:,2,3) + De_E(:,:,3,2));
    De_S(:,:,3,2) = De_S(:,:,2,3);
else % anisotropic material
    for yi = 1:size(A, 1)
        for xi = 1:size(A, 2)
            e_voight = [De_E(yi, xi, 1, 1); De_E(yi, xi, 2, 2); De_E(yi, xi, 3, 3); ...
                De_E(yi, xi, 1, 2) + De_E(yi, xi, 2, 1); ...
                De_E(yi, xi, 1, 3) + De_E(yi, xi, 3, 1); ...
                De_E(yi, xi, 2, 3) + De_E(yi, xi, 3, 2)];
            s_voight = Maps.Stiffness * e_voight;
            De_S(yi, xi, :, :) = [s_voight(1), s_voight(4), s_voight(5); ...
                s_voight(4), s_voight(2), s_voight(6); ...
                s_voight(5), s_voight(6), s_voight(3)];
            De_E(yi, xi, :, :) = [e_voight(1), e_voight(4)/2, e_voight(5)/2; ...
                e_voight(4)/2, e_voight(2), e_voight(6)/2; ...
                e_voight(5)/2, e_voight(6)/2, e_voight(3)];
        end
    end
end
end

%%
function [ alldata,dataum ] = reshapeDefromationGradient( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
x  = raw_data(:,1);
y  = raw_data(:,2);
z  = raw_data(:,3);
xVec = unique(x);
yVec = unique(y);
zVec = unique(z);
du11 = raw_data(:,4);
du12 = raw_data(:,5);
du13 = raw_data(:,6);
du21 = raw_data(:,7);
du22 = raw_data(:,8);
du23 = raw_data(:,9);
du31 = raw_data(:,10);
du32 = raw_data(:,11);
du33 = raw_data(:,12);

[dataum.X,dataum.Y,dataum.Z] = meshgrid(xVec,yVec,zVec);
[nRows, nCols , nDep] = size(dataum.X);
dataum.du11 = zeros(nRows, nCols, nDep); %Initialise
dataum.du12 = zeros(nRows, nCols, nDep); %Initialise
dataum.du13 = zeros(nRows, nCols, nDep); %Initialise
dataum.du21 = zeros(nRows, nCols, nDep); %Initialise
dataum.du22 = zeros(nRows, nCols, nDep); %Initialise
dataum.du23 = zeros(nRows, nCols, nDep); %Initialise
dataum.du31 = zeros(nRows, nCols, nDep); %Initialise
dataum.du32 = zeros(nRows, nCols, nDep); %Initialise
dataum.du33 = zeros(nRows, nCols, nDep); %Initialise

for iRow = 1:nRows % loop rows
    for iCol = 1:nCols % loop cols
        for iDep = 1:nDep
            xt = dataum.X(iRow,iCol,iDep);
            yt = dataum.Y(iRow,iCol,iDep);
            zt = dataum.Z(iRow,iCol,iDep);
            idx = find(x==xt & y==yt & z==zt);
            if ~isempty(idx)
                du11t = du11(idx(1));
                du12t = du12(idx(1));
                du13t = du13(idx(1));
                du21t = du21(idx(1));
                du22t = du22(idx(1));
                du23t = du23(idx(1));
                du31t = du31(idx(1));
                du32t = du32(idx(1));
                du33t = du33(idx(1));
                dataum.du11(iRow,iCol,iDep) = du11t;
                dataum.du12(iRow,iCol,iDep) = du12t;
                dataum.du13(iRow,iCol,iDep) = du13t;
                dataum.du21(iRow,iCol,iDep) = du21t;
                dataum.du22(iRow,iCol,iDep) = du22t;
                dataum.du23(iRow,iCol,iDep) = du23t;
                dataum.du31(iRow,iCol,iDep) = du31t;
                dataum.du32(iRow,iCol,iDep) = du32t;
                dataum.du33(iRow,iCol,iDep) = du33t;
            end
        end
    end
end
alldata = [dataum.X(:)      dataum.Y(:)      dataum.Z(:)      dataum.du11(:) ...
    dataum.du12(:)   dataum.du13(:)   dataum.du21(:)   dataum.du22(:)   ...
    dataum.du23(:)   dataum.du31(:)   dataum.du32(:)   dataum.du33(:)];

end
%%
function [E,v,G,Co] = effectiveE_nu(C)
% this function caclulates effective Youn modulus and Possion ratio for a
% ansitropic material based on this paper
% Reference: https://doi.org/10.3390/cryst8080307
BV = (C(1,1)+2*C(1,2))/3;               % Voigt bulk modulus
GV = (C(1,1)-C(1,2)+3*C(4,4))/5;        % Voigt shear modulus

S = C^-1;
BR = 1/(3*S(1,1)+6*S(1,2));             % Reuss bulk modulus
GR = 5/(4*S(1,1)-4*S(1,2)+3*S(4,4));    % Reuss shear modulus

B = (BR+BV)/2;                          % Hills average bulk modulus
G = (GR+GV)/2;                          % Hills average shear modulus
E = 9*B*G/(3*B+G);                      % Youngs modulus (E)
v = (3*B-E)/(6*B);                      % Poissons ratio
Co = [];

% K = (C(1,1)+C(2,2)+C(3,3)+2*(C(1,2)+C(2,3)+C(1,2)))/9; % istropic shear Modulus
% Gv = (C(1,1)+C(2,2)+C(3,3)-(C(1,2)+C(2,3)+C(1,2))+2*(C(4,4)+C(5,5)+C(6,6)))/15; % Bulk Modulus

%% Paper: What is the Youngs Modulus of Silicon?
Cc =C^-1;
Co.Ex = 1/Cc(1,1);
Co.Ey = 1/Cc(2,2);
Co.Ez = 1/Cc(3,3);
Co.Gxy = 1/Cc(4,4);
Co.Gxz = 1/Cc(5,5);
Co.Gyz = 1/Cc(6,6);
Co.vxy = -Co.Ey*Cc(1,2);
Co.vxz = -Co.Ez*Cc(1,3);
Co.vyz = -Co.Ez*Cc(2,3);

Co.C = [ 1/Co.Ex       -Co.vxy/Co.Ey   -Co.vxz/Co.Ez 0   0   0
    -Co.vxy/Co.Ex    1/Co.Ey      -Co.vyz/Co.Ez 0   0   0
    -Co.vxz/Co.Ex   -Co.vyz/Co.Ey  1/Co.Ez     0   0   0
    0          0           0 	1/Co.Gyz 0   0
    0          0           0   0   1/Co.Gxz 0
    0          0           0   0   0   1/Co.Gxy];
Co.C = Co.C^-1;

%% a different approach as sometime the first approach sometimes
% delivers minus results!
if G<0 || E<0 || v<0 || v > 0.5
    if v > 0.5 || v<0 % for metals
        v = abs((3*BV-E)/(6*BV));
    end
end

end
%%
function [Crop] = CroppingEqually(Maps)
close all;                  fig=subplot(1,1,1);
imagesc(Maps.X(1,:),Maps.Y(:,1),Maps.du11);
axis image;                 set(gca,'Ydir','normal');   %axis off;
colorbar;   %colormap jet;
set(gcf,'position',[30 50 1300 950])
xlabel('X [Raw Data Units]');          ylabel('Y [Raw Data Units]');
title('select the tip')
[lineX,lineY] = ginput(1);
title('Select Area to Crop');
[Xcrop,Ycrop] = ginput(2);
Xcrop = [min(Xcrop);max(Xcrop)];
Ycrop = [min(Ycrop);max(Ycrop)];
xLin          = Maps.X(1,:);                     yLin         = Maps.Y(:,1);
[~, Xcrop(1)] = min(abs(xLin-Xcrop(1)));   Xcrop(1) = xLin(Xcrop(1));
[~, Xcrop(2)] = min(abs(xLin-Xcrop(2)));   Xcrop(2) = xLin(Xcrop(2));
[~, Ycrop(1)] = min(abs(yLin-Ycrop(1)));   Ycrop(1) = yLin(Ycrop(1));
[~, Ycrop(2)] = min(abs(yLin-Ycrop(2)));   Ycrop(2) = yLin(Ycrop(2));
if abs(mean(lineY)-Ycrop(1)) ~= abs(mean(lineY)-Ycrop(2))
    addi  = (abs(lineY-Ycrop(1))+abs(lineY-Ycrop(2)))/2;
    Ycrop = [lineY-addi, lineY+addi];
end
Xcrop(1) = 2*lineX-Xcrop(2);
Dis = (abs(Ycrop(2) - Ycrop(1))-abs(Xcrop(2) - Xcrop(1)))/2;
Xcrop = [Xcrop(1)-Dis Xcrop(2)+Dis];
hold on
plot([Xcrop(1) Xcrop(2) Xcrop(2) Xcrop(1) Xcrop(1)],...
    [Ycrop(1) Ycrop(1) Ycrop(2) Ycrop(2) Ycrop(1)],'color','k')
plot(lineX,lineY,'pw')
hold off

[~, Xcrop(1)] = min(abs(xLin-Xcrop(1)));
[~, Xcrop(2)] = min(abs(xLin-Xcrop(2)));
[~, Ycrop(1)] = min(abs(yLin-Ycrop(1)));
[~, Ycrop(2)] = min(abs(yLin-Ycrop(2)));

for iV=1:3
    for iO=1:3
        eval(sprintf('Crop.du%d%d = Maps.du%d%d(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));',...
            iV,iO,iV,iO));
    end
end

%% XY, steps and stifness
Maps.Z = zeros(size(Maps.X));
Crop.X   = Maps.X(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.Y   = Maps.Y(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.Z   = Maps.Z(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.X   = Crop.X - min(min(Crop.X));  	Crop.Y   = Crop.Y - min(min(Crop.Y));
if (Crop.X(1) - Crop.X(end))>0;         Crop.X   = flip(Crop.X,2);         end
if (Crop.Y(1) - Crop.Y(end))>0;         Crop.Y   = flip(Crop.Y,1);         end

end

%%
function plot_DecomposeddU(du_dx,Maps)
figure;
iV=0;   Mo = {'I','II','III'}; KK = {'x','y','z'};
for ii=1:3
    for ij=1:3
        eval(sprintf('pD = Maps.du%d%d;',ii,ij));
        if ~sum(pD(:)) == 0
            iV= iV+1;
            s{iV}=subplot(9,3,iV);
            pcolor(Maps.X,Maps.Y,pD); clear pD
            title(['\nablau_{' KK{ii} KK{ij} '}'],'fontsize',12);   shading interp;
            axis image; axis off; colormap jet;         box off;
            c  =colorbar;	cU(iV,:) = c.Limits;         colorbar off;
        end
    end
end
for iO =1:3
    for ii=1:3
        for ij=1:3
            pD = squeeze(du_dx(:,:,ii,ij,iO));
            if ~sum(pD(:)) == 0
                iV= iV+1;
                s{iV}=subplot(9,3,iV);
                pcolor(Maps.X,Maps.Y,pD); clear pD
                title(['\nablau^{' Mo{iO} '}_{' KK{ii} KK{ij} '}'],'fontsize',12);
                shading interp;
                axis image; axis off; colormap jet;         box off;
                c  =colorbar;	cU(iV,:) = c.Limits;         colorbar off;
            end
        end
    end
end
addScale([9 3 iV],[Maps.X(:) Maps.Y(:)]);

cbax  = axes('visible', 'off');             cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax,'location','westoutside','position',[0.9011 0.1211 0.0121 0.7533]);
h.Label.String = '\nablau';
h.Label.FontSize = 30;
for iO = 1:length(s)
    set(s{iO},"clim",caxis);
end
%}
set(gcf,'position',[348 59 1396 932]);
end

%%
function plot_DecomposedStrain(uXXd,uYYd,uZZd,uXYd,uXZd,uYZd,Maps)
figure;
s1=subplot(3,3,1);  	pcolor(Maps.X,Maps.Y,Maps.E11);
title([char(949) '_{xx}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off;
c  =colorbar;	cU(1,:) = c.Limits;         colorbar off;
s2=subplot(3,3,2);  	pcolor(Maps.X,Maps.Y,Maps.E12);
title([char(949) '_{xy}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;         colorbar off;
s3=subplot(3,3,3);  	pcolor(Maps.X,Maps.Y,Maps.E31);
title([char(949) '_{xz}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;         colorbar off;
s5=subplot(3,3,5);  	pcolor(Maps.X,Maps.Y,Maps.E22);
title([char(949) '_{yy}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;         colorbar off;
s6=subplot(3,3,6);  	pcolor(Maps.X,Maps.Y,Maps.E32);
title([char(949) '_{yz}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;         colorbar off;
s9=subplot(3,3,9);  	pcolor(Maps.X,Maps.Y,Maps.E33);
title([char(949) '_{zz}'],'fontsize',19);   shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;         colorbar off;
addScale([3 3 9],[Maps.X(:) Maps.Y(:)]);

EId   = sqrt(0.5*((uXXd(:,:,1)-uYYd(:,:,1)).^2+(uYYd(:,:,1)-uZZd(:,:,1)).^2+...
    (uZZd(:,:,1)-uXXd(:,:,1)).^2+ ...
    (uXYd(:,:,1).^2+uYZd(:,:,1).^2+uXZd(:,:,1).^2).*6));
EIId  = sqrt(0.5*((uXXd(:,:,2)-uYYd(:,:,2)).^2+(uYYd(:,:,2)-uZZd(:,:,2)).^2+...
    (uZZd(:,:,2)-uXXd(:,:,2)).^2+ ...
    (uXYd(:,:,2).^2+uYZd(:,:,2).^2+uXZd(:,:,2).^2).*6));
EIIId = sqrt(0.5*((uXXd(:,:,3)-uYYd(:,:,3)).^2+(uYYd(:,:,3)-uZZd(:,:,3)).^2+...
    (uZZd(:,:,3)-uXXd(:,:,3)).^2+ ...
    (uXYd(:,:,3).^2+uYZd(:,:,3).^2+uXZd(:,:,3).^2).*6));

s4=subplot(3,3,4);  	pcolor(Maps.X,Maps.Y,EId);
title([char(949) '^{I}_M'],'fontsize',19);  shading interp;
axis image; axis off;  box off;             colormap jet;
c  =colorbar;	cU(7,:) = c.Limits;         colorbar off;
s7=subplot(3,3,7);  	pcolor(Maps.X,Maps.Y,EIId);
title([char(949) '^{II}_M'],'fontsize',19); shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(8,:) = c.Limits;         colorbar off;
s8=subplot(3,3,8);  	pcolor(Maps.X,Maps.Y,EIIId);
title([char(949) '^{III}_M'],'fontsize',19);shading interp;
axis image; axis off; colormap jet;         box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(9,:) = c.Limits;         colorbar off;
%
cbax  = axes('visible', 'off');             cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax,'location','westoutside','position',[0.9011 0.1211 0.0121 0.7533]);
h.Label.String = [char(949)];
h.Label.FontSize = 30;
set([s1 s2 s3 s4 s5 s6 s7 s8 s9],"clim",caxis);
%}
set(gcf,'position',[348 59 1396 932]);
end

%%
function plot_DecomposedStess(uXXd,uYYd,uZZd,uXYd,uXZd,uYZd,Maps,Saf)
figure;
s1=subplot(3,3,1);  	pcolor(Maps.X,Maps.Y,Maps.S11*Saf*1e-9);
title('\sigma_{xx}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off;
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off;
s2=subplot(3,3,2);  	pcolor(Maps.X,Maps.Y,Maps.S12*Saf*1e-9);
title('\sigma_{xy}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	pcolor(Maps.X,Maps.Y,Maps.S13*Saf*1e-9);
title('\sigma_{xz}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	pcolor(Maps.X,Maps.Y,Maps.S22*Saf*1e-9);
title('\sigma_{yy}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	pcolor(Maps.X,Maps.Y,Maps.S23*Saf*1e-9);
title('\sigma_{yz}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;     colorbar off;
s9=subplot(3,3,9);  	pcolor(Maps.X,Maps.Y,Maps.S33*Saf*1e-9);
title('\sigma_{zz}','fontsize',19);     shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
addScale([3 3 9],[Maps.X(:) Maps.Y(:)]);

SId   = sqrt(0.5*((uXXd(:,:,1)-uYYd(:,:,1)).^2+(uYYd(:,:,1)-uZZd(:,:,1)).^2+...
    (uZZd(:,:,1)-uXXd(:,:,1)).^2+ ...
    (uXYd(:,:,1).^2+uYZd(:,:,1).^2+uXZd(:,:,1).^2).*6));
SIId  = sqrt(0.5*((uXXd(:,:,2)-uYYd(:,:,2)).^2+(uYYd(:,:,2)-uZZd(:,:,2)).^2+...
    (uZZd(:,:,2)-uXXd(:,:,2)).^2+ ...
    (uXYd(:,:,2).^2+uYZd(:,:,2).^2+uXZd(:,:,2).^2).*6));
SIIId = sqrt(0.5*((uXXd(:,:,3)-uYYd(:,:,3)).^2+(uYYd(:,:,3)-uZZd(:,:,3)).^2+...
    (uZZd(:,:,3)-uXXd(:,:,3)).^2+ ...
    (uXYd(:,:,3).^2+uYZd(:,:,3).^2+uXZd(:,:,3).^2).*6));

s4=subplot(3,3,4);  	pcolor(Maps.X,Maps.Y,SId*1e-9);
title('\sigma^{I}_M','fontsize',19);    shading interp;
axis image; axis off;  box off;         colormap jet;
c  =colorbar;	cU(7,:) = c.Limits;     colorbar off;
s7=subplot(3,3,7);  	pcolor(Maps.X,Maps.Y,SIId*1e-9);
title('\sigma^{II}_M','fontsize',19);   shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(8,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	pcolor(Maps.X,Maps.Y,SIIId*1e-9);
title('\sigma^{III}_M','fontsize',19);  shading interp;
axis image; axis off; colormap jet;     box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(9,:) = c.Limits;     colorbar off;
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax,'location','westoutside','position',[0.9011 0.1211 0.0121 0.7533]);
h.Label.String = '\sigma [GPa]';
h.Label.FontSize = 30;
set([s1 s2 s3 s4 s5 s6 s7 s8 s9],"clim",caxis);
%}
set(gcf,'position',[348 59 1396 932]);
end


%%
function plot_JKIII(KI,KII,KIII,J,stepsize,input_unit)
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',14)
Contour = (1:length(J.Raw))*stepsize;
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;    hold on;
plot(Contour,KI.Raw,'k--o','MarkerEdgeColor','k','LineWidth',4);
plot(Contour,KII.Raw,'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
plot(Contour,KIII.Raw,'k--d','MarkerEdgeColor','k','LineWidth',4');
ylabel('K (MPa m^{0.5})'); hold off
Kd = [KI.Raw(:); KII.Raw(:); KIII.Raw(:)];
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
yyaxis right;
plot(Contour,J.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
ylabel('J (J/m^2)');        ylim([0 max(J.Raw)+min(J.Raw)/4]);
legend(['K_{I} = '     num2str(KI.true)   ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
    ['K_{II} = '       num2str(KII.true)  ' ± ' num2str(KII.div) ' MPa\surdm' ],...
    ['K_{III} = '      num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
    ['J_{integral} = ' num2str(J.true)    ' ± ' num2str(J.div)   ' J/m^2'],...
    'location','northoutside','box','off');
set(gcf,'position',[60,-70,750,1000]);grid on;  box off;
ax1 = gca;  axPos = ax1.Position;
% Change the position of ax1 to make room for extra axes
% format is [left bottom width height], so moving up and making shorter here...
ax1.Position = axPos + [0 0.2 0 -0.15];
% Exactly the same as for plots (above), axes LineWidth can be changed inline or after
ax1.LineWidth = 1;
% Add two more axes objects, with small multiplier for height, and offset for bottom
ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 0.08 0 0], 'color', 'none', 'linewidth', 1);
% You can change the limits of the new axes using XLim
ax2.XLim = [0 length(Contour)+1];     ax1.XLim = [0 max(Contour)+stepsize];
% You can label the axes using XLabel.String
if strcmpi(input_unit,'um')
    input_unit = '\mum';
end
ax1.XLabel.String = ['Contour Distance [' input_unit ']'];
ax2.XLabel.String = 'Contour Number';
end

%%
function addScale(No,alldata)
% funciton to add a measurment line to the graph, you can either input
% number of the suplots or the exact suplot where you want to add the scale
% bar
if length(No) == 1      && No == 2
    subplot(1,2,1); AddScalePar(alldata(:,1),alldata(:,2))
    subplot(1,2,2); AddScalePar(alldata(:,1),alldata(:,2))
elseif length(No) == 1  && No == 3
    subplot(1,3,1); AddScalePar(alldata(:,1),alldata(:,2))
    subplot(1,3,2); AddScalePar(alldata(:,1),alldata(:,2))
    subplot(1,3,3); AddScalePar(alldata(:,1),alldata(:,2))
elseif length(No) == 1	&& No == 1
    AddScalePar(alldata(:,1),alldata(:,2));
elseif length(No) == 3
    subplot(No(1),No(2),No(3)); AddScalePar(alldata(:,1),alldata(:,2));
end
end

function AddScalePar(X,Y)
X = unique(X);         Y = unique(Y);
hold on; line([X(ceil(end-length(X)*0.05)),X(end-ceil(length(X)*0.15))],...
    [Y(ceil(length(Y)*0.05)),Y(ceil(length(Y)*0.05))],'Color','k','LineWidth',5)
ht = text(double(X(end-ceil(length(X)*0.24))),double(Y(ceil(length(Y)*0.14)))...
    ,[num2str(round(abs(X(ceil(length(X)*0.05))-X(ceil(length(X)*0.15))),1)) '\mum']);
set(ht,'FontSize',16,'FontWeight','Bold')
set(gca,'Visible','off')
hold off
end

%%
function plot_JM(M, J, stepsize, input_unit,LorM)
Md = M.Raw(:);
set(0, 'defaultAxesFontSize', 22);
set(0, 'DefaultLineMarkerSize', 14);
Contour = (1:length(J.vectorial)) * stepsize;

% Create figure and set color order
fig = figure;
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);

% Plot M on the left y-axis
yyaxis left;
hold on;
% plot(Contour, M.Raw, 'k--d', 'MarkerEdgeColor', 'k', 'LineWidth', 4);  % Plot M.MRaw
plot(Contour, M.Raw(1,:), 'k--o', 'LineWidth', 4, 'MarkerFaceColor','k');  % Plot M.MRaw
plot(Contour, M.Raw(2,:), 'k--s', 'LineWidth', 4, 'MarkerFaceColor','k');  % Plot M.MRaw
ylabel([LorM ' (J/m)']);
Kd = M.Raw(:);
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
hold off;

% Plot J on the right y-axis
yyaxis right; hold on
plot(Contour, J.vectorial(1,:), 'r-->', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');  % Plot J.JRaw
plot(Contour, J.vectorial(2,:), 'r--<', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');  % Plot J.JRaw

ylabel('J (J/m^2)');
Kd = J.vectorial(:);
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
hold off

% Update legend to only include M and J
% legend(['M_{integral} = ' num2str(M.true) ' ± ' num2str(M.div) ' J/m'], ... %N/m
legend([LorM '_{1} = ' num2str(M.true(1)) ' ± ' num2str(M.div(1)) ' J/m'], ...
    [LorM '_{2} = ' num2str(M.true(2)) ' ± ' num2str(M.div(2)) ' J/m'], ...
    ['J_{1} = ' num2str(J.vectorial_true(1)) ' ± ' num2str(J.vectorial_div(1)) ' J/m^2'], ... 
    ['J_{2} = ' num2str(J.vectorial_true(2)) ' ± ' num2str(J.vectorial_div(2)) ' J/m^2'], ...%N
    'location', 'northoutside', 'box', 'off');

% Set figure size and formatting
set(gcf, 'position', [800, -70, 750, 1000]);
grid on;
box off;

% Adjust axes position
ax1 = gca;
axPos = ax1.Position;
ax1.Position = axPos + [0 0.2 0 -0.15];
ax1.LineWidth = 1;

% Add second x-axis for contour number
ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 0.08 0 0], ...
    'color', 'none', 'linewidth', 1);
ax2.XLim = [0 length(Contour) + 1];
ax1.XLim = [0 max(Contour) + stepsize];

% Label the axes
if strcmpi(input_unit, 'um')
    input_unit = '\mum';
end
ax1.XLabel.String = ['Contour Distance [' input_unit ']'];
ax2.XLabel.String = 'Contour Number';
end

%%
function [ alldata,dataum ] = reshapeData( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
% try
%     if size(raw_data,2) == 4
%         dataum.X1  = reshape(raw_data(:,1),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Y1  = reshape(raw_data(:,2),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Ux  = reshape(raw_data(:,3),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         dataum.Uy  = reshape(raw_data(:,4),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))));
%         alldata = [dataum.X1(:) dataum.Y1(:) dataum.Ux(:) dataum.Uy(:)];
%     else
%         dataum.X1  = reshape(raw_data(:,1),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Y1  = reshape(raw_data(:,2),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Z1  = reshape(raw_data(:,3),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Ux  = reshape(raw_data(:,4),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Uy  = reshape(raw_data(:,5),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         dataum.Uz  = reshape(raw_data(:,6),length(unique(raw_data(:,1))),length(unique(raw_data(:,2))),length(unique(raw_data(:,3))));
%         alldata = [dataum.X1(:) dataum.Y1(:) dataum.Z1(:) dataum.Ux(:) dataum.Uy(:) dataum.Uz(:)];
%     end
% catch
    
    x  = raw_data(:,1);
    y  = raw_data(:,2);
    ux = raw_data(:,3);
    uy = raw_data(:,4);
    
    xVec = unique(x);
    yVec = unique(y);
    
    % nDataPoints = length(x);
    % if length(xVec) <	length(raw_data(:,1))*0.5
    
    %Define grid
    [xMap,yMap] = meshgrid(xVec,yVec);
    [nRows, nCols] = size(xMap);
    
    % nGridPoints = length(xMap(:));
    
    uxMap = NaN(nRows, nCols); %Initialise
    uyMap = NaN(nRows, nCols); %Initialise
    
    if size(raw_data,2) == 6
        z  = raw_data(:,3);
        zVec = unique(z);
        ux = raw_data(:,4);
        uy = raw_data(:,5);
        uz = raw_data(:,6);
        [xMap,yMap,zMap] = meshgrid(xVec,yVec,zVec);
        [nRows, nCols , nDep] = size(xMap);
        uxMap = NaN(nRows, nCols, nDep); %Initialise
        uyMap = NaN(nRows, nCols, nDep); %Initialise
        uzMap = NaN(nRows, nCols, nDep); %Initialise
    end
    
    for iRow = 1:nRows % loop rows
        for iCol = 1:nCols % loop cols
            
            if size(raw_data,2) == 6 %% 3D
                for iDep = 1:nDep
                    xt = xMap(iRow,iCol,iDep);
                    yt = yMap(iRow,iCol,iDep);
                    zt = zMap(iRow,iCol,iDep);
                    idx = find(x==xt & y==yt & z==zt);
                    if ~isempty(idx)
                        uxt = ux(idx(1));
                        uyt = uy(idx(1));
                        uzt = uz(idx(1));
                        uxMap(iRow,iCol,iDep) = uxt;
                        uyMap(iRow,iCol,iDep) = uyt;
                        uzMap(iRow,iCol,iDep) = uzt;
                    end
                end
                
            else %% 2D
                xt = xMap(iRow,iCol);
                yt = yMap(iRow,iCol);
                idx = find(and(x==xt,y==yt)); %find linear index of point corresponding to xt,yt;
                if ~isempty(idx)
                    uxt = ux(idx(1));
                    uyt = uy(idx(1));
                    uxMap(iRow,iCol) = uxt;
                    uyMap(iRow,iCol) = uyt;
                end
            end
        end
    end
    
    dataum.X1 = xMap;
    dataum.Y1 = yMap;
    % dataum.Uy = uyMap;
    % dataum.Ux = uxMap;
    % threshold = 0.95;
    % [ uxMap ] = dispFieldSmoothing( uxMap, threshold );
    dataum.Ux = uxMap;
    % [ uyMap ] = dispFieldSmoothing( uyMap, threshold );
    dataum.Uy = uyMap;

alldata = [dataum.X1(:) dataum.Y1(:) dataum.Ux(:) dataum.Uy(:)];
if size(raw_data,2) == 6
    dataum.Z1 = zMap;
    dataum.Uz = uzMap;
    alldata = [dataum.X1(:) dataum.Y1(:) dataum.Z1(:) dataum.Ux(:) dataum.Uy(:) dataum.Uz(:)];
end
% end
%{
else
	disp('the data is highly non-uniform, I will now re-arrange the data');
%     scatter3(alldata(:,1), alldata(:,2), alldata(:,3),[],alldata(:,3)); view([0 90])
	Fx = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,3),'natural','nearest');
	Fy = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,4),'natural','nearest');
	X = linspace(min(raw_data(:,1)),max(raw_data(:,1)),300);
	Y = min(raw_data(:,2)):abs(X(2)-X(1)):max(raw_data(:,2));
	[dataum.X1,dataum.Y1] = meshgrid(X,Y);
	Ux = Fx(dataum.X1(:),dataum.Y1(:));
	Uy = Fy(dataum.X1(:),dataum.Y1(:));
    dataum.Ux = reshape(Ux,length(Y),length(X));
    dataum.Uy = reshape(Uy,length(Y),length(X));
	% scatter3(x(:),y(:),Ux,[],Ux); view([0 90])
end
%}
end

function [cx,cy,cz]=crackgradient(c,dx)
c=squeeze(c);
[row,~]=size(c);
midr=floor(row/2);
ctop=c(1:midr,:);
cbot=c(midr+1:end,:);
[cxtop,cytop]=gradient(ctop,dx);
[cxbot,cybot]=gradient(cbot,dx);
cx=[cxtop;cxbot];
cy=[cytop;cybot];
cz=zeros(size(cx));
end
