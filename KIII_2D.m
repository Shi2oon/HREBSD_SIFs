function [J,KI,KII,KIII] = KIII_2D(Maps,MatProp,xEBSD)
close all; clc

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

% Another option is to input the strain components as a vector matrix with
% 9 columns the first three columns are the x, y and z coordinate in meters. z
% cooridnate can be a zero column. the 4th to the 9th column are the strain
% components arranged as
% Maps = [X(:) Y(:) Z(:) E11(:) E12(:) E13(:) E21(:) E22(:) E23(:) E31(:) E32(:) E33(:)]; 

% if the map is a 2D strain map then zero all out of the plane components

% the material paramters as MatProp.E for Young's Modulus and MatProp.nu
% for Possions ratio. or as a stifness matrix all in Pa
% xEBSD will assume and solve for  plane strain conditions with sigma33 == 0 + no
% volumetric change.

%
if size(Maps,2) > 1
    alldata = Maps; clear Maps
    [~,Maps]=reshapeStrain(alldata);
    Maps.stepsize = unique(round(diff(unique(Maps.Y(:))),3));
    Maps.units.xy = 'm';
    Maps.units.St = 'Pa';
    if size(MatProp,1) == 6
        Maps.Stiffness = MatProp;
    else
        Maps.E  = MatProp.E;
        Maps.nu = MatProp.nu;
    end
end

%% prepare Data
imagesc(Maps.E11);axis tight; axis image; axis off
set(gcf,'position',[737 287 955 709]);
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'N';         % Use the TeX interpreter to format the question
quest            = 'Do you want to Crop and Centre the Crack tip';
answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if strcmpi(answer,'Y') % crop data
    [Crop] = CroppingEqually(Maps);
    Maps.X   = Crop.X;      Maps.Y   = Crop.Y;      Maps.Z   = Crop.Z;
    Maps.E11 = Crop.E11;    Maps.E12 = Crop.E12;    Maps.E13 = Crop.E13;
    Maps.E21 = Crop.E21;    Maps.E22 = Crop.E22;    Maps.E23 = Crop.E23;
    Maps.E31 = Crop.E31;    Maps.E32 = Crop.E32;    Maps.E33 = Crop.E33;
end
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'L';         % Use the TeX interpreter to format the question
quest            = 'Is the crack on your left or right ?';
answer           = questdlg(quest,'Boundary Condition','L','R', opts);
if strcmpi(answer,'R') % crop data
    Maps.E11 = flip(flip(Maps.E11,1),2);    Maps.E12 = flip(flip(Maps.E12,1),2);
    Maps.E13 = flip(flip(Maps.E13,1),2);
    Maps.E21 = flip(flip(Maps.E21,1),2);    Maps.E22 = flip(flip(Maps.E22,1),2);
    Maps.E23 = flip(flip(Maps.E23,1),2);
    Maps.E31 = flip(flip(Maps.E31,1),2);    Maps.E32 = flip(flip(Maps.E32,1),2);
    Maps.E33 = flip(flip(Maps.E33,1),2);
end

% close
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
Maps.stepsize = Maps.stepsize*saf;
Maps.units.St = 'Pa';        Maps.units.xy = 'm';

%% Input
DataSize = [size(Maps.E11),1];

%% JMAN approach (without FEM) - Standard J-integral.
% An approach to calculate the J -integral by digital image correlation
% displacement field measurement, FFEMS (2012), 35, ?971-984
% Displacement gradient
uXX = Maps.E11;     uXY = Maps.E12;     uXZ = Maps.E13;
uYX = Maps.E21;     uYY = Maps.E22;     uYZ = Maps.E23;
uZX = Maps.E31;     uZY = Maps.E32;     uZZ = Maps.E33;
% Infitisimal strain tensor
E11 = uXX;
E22 = uYY;
E33 = uZZ;
E12 = 1/2*(uXY+uYX);
E13 = 1/2*(uXZ+uZX);
E23 = 1/2*(uYZ+uZY);

if isfield(Maps,'Stiffness')
    Maps.Stiffness = Maps.Stiffness*Saf;
    [Maps.E,Maps.nu,Maps.G,Maps.Co] = effectiveE_nu(Maps.Stiffness); % in Pa
    for iR = 1:size(Maps.E11,1)
        for iC = 1:size(Maps.E11,2)
            strain(iR,iC,:) = [uXX(iR,iC);   uYY(iR,iC);    uZZ(iR,iC);...
                2*E12(iR,iC); 2*E13(iR,iC);  2*E23(iR,iC)];
            if exist('xEBSD','var')
                %solve the (HR-EBSD) boundary condition
                K1=strain(iR,iC,1)-strain(iR,iC,3);
                K2=strain(iR,iC,2)-strain(iR,iC,3);
                K3=strain(iR,iC,4)*Maps.Stiffness(4,3)+...
                    strain(iR,iC,5)*Maps.Stiffness(5,3)+...
                    strain(iR,iC,6)*Maps.Stiffness(6,3);
                
                e33n=-(K1*Maps.Stiffness(1,3)+K2*Maps.Stiffness(2,3)+K3)/...
                    (Maps.Stiffness(1,3)+Maps.Stiffness(2,3)+Maps.Stiffness(3,3));
                e11n=K1+e33n;
                e22n=K2+e33n;
                %new strain/stress vector
                co = strain(iR,iC,4:6);
                strain(iR,iC,:)=[e11n;e22n;e33n;co(:)];
            end
            stress(iR,iC,:) = Maps.Stiffness*squeeze(strain(iR,iC,:));
            sXX(iR,iC) = stress(iR,iC,1);
            sYY(iR,iC) = stress(iR,iC,2);
            sZZ(iR,iC) = stress(iR,iC,3);
            sXY(iR,iC) = stress(iR,iC,4);
            sXZ(iR,iC) = stress(iR,iC,5);
            sYZ(iR,iC) = stress(iR,iC,6);
        end
    end
    E11 = squeeze(strain(:,:,1));
    E22 = squeeze(strain(:,:,2));
    E33 = squeeze(strain(:,:,3));
    E12 = 0.5*squeeze(strain(:,:,4));
    E13 = 0.5*squeeze(strain(:,:,5));
    E23 = 0.5*squeeze(strain(:,:,6));
else    % Chauchy stress tensor, assming linear-elastic, isotropic material for
    % plane stress conditions
    Maps.E = Maps.E*Saf;
    Maps.G = Maps.E/(2*(1 + Maps.nu));
    sXX = Maps.E/(1-Maps.nu^2)*(E11+Maps.nu*(E22+E33));
    sYY = Maps.E/(1-Maps.nu^2)*(E22+Maps.nu*(E11+E33));
    sZZ = Maps.E/(1-Maps.nu^2)*(E33+Maps.nu*(E11+E22));
    sXY = 2*Maps.G*E12;
    sXZ = 2*Maps.G*E13;
    sYZ = 2*Maps.G*E23;
end
clear stress strain
% Elastic strain energy
W = 0.5*(E11.*sXX+E22.*sYY+E33.*sZZ+2*E12.*sXY+2*E13.*sXZ+2*E23.*sYZ);
% Generate q field
celw = 1;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(Maps.stepsize*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*Maps.stepsize^2;
JA = ((sXX.*uXX+sXY.*uYX+sXZ.*uZX-W).*dQdX+(sYY.*uYX+sXY.*uXX+sYZ.*uZX).*dQdY).*dA;
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2));
    J.JIII(ii) = sum(sum(JA.*areaID,'omitnan'),'omitnan');
end
J.Raw = abs(J.JIII);

%% Decomposition method.
% Mode I
uXXd(:,:,1) = 0.5*(Maps.E11+flipud(Maps.E11));
uXYd(:,:,1) = 0.5*(Maps.E12-flipud(Maps.E12));
uXZd(:,:,1) = 0.5*(Maps.E13+flipud(Maps.E13));

uYXd(:,:,1) = 0.5*(Maps.E21-flipud(Maps.E21));
uYYd(:,:,1) = 0.5*(Maps.E22+flipud(Maps.E22));
uYZd(:,:,1) = 0.5*(Maps.E23-flipud(Maps.E23));

uZXd(:,:,1) = 0.5*(Maps.E31+flipud(Maps.E31));
uZYd(:,:,1) = 0.5*(Maps.E32-flipud(Maps.E32));
uZZd(:,:,1) = 0.5*(Maps.E33+flipud(Maps.E33));

% Mode II
uXXd(:,:,2) = 0.5*(Maps.E11-flipud(Maps.E11));
uXYd(:,:,2) = 0.5*(Maps.E12+flipud(Maps.E12));
uXZd(:,:,2) = zeros(DataSize);

uYXd(:,:,2) = 0.5*(Maps.E21+flipud(Maps.E21));
uYYd(:,:,2) = 0.5*(Maps.E22-flipud(Maps.E22));
uYZd(:,:,2) = zeros(DataSize);

uZXd(:,:,2) = zeros(DataSize);
uZYd(:,:,2) = zeros(DataSize);
uZZd(:,:,2) = 0.5*(Maps.E33-flipud(Maps.E33));

% Mode III
uXXd(:,:,3) = zeros(DataSize);
uXYd(:,:,3) = zeros(DataSize);
uXZd(:,:,3) = 0.5*(Maps.E13-flipud(Maps.E13));

uYXd(:,:,3) = zeros(DataSize);
uYYd(:,:,3) = zeros(DataSize);
uYZd(:,:,3) = 0.5*(Maps.E23+flipud(Maps.E23));

uZXd(:,:,3) = 0.5*(Maps.E31-flipud(Maps.E31));
uZYd(:,:,3) = 0.5*(Maps.E32+flipud(Maps.E32));
uZZd(:,:,3) = 0.5*(Maps.E33-flipud(Maps.E33));

% Infitisimal strain tensor
E11d = uXXd;
E22d = uYYd;
E33d = uZZd;
E12d = 1/2*(uXYd+uYXd);
E13d = 1/2*(uXZd+uZXd);
E23d = 1/2*(uYZd+uZYd);

if isfield(Maps,'Stiffness')
    for iR = 1:size(Maps.E11,1)
        for iC = 1:size(Maps.E11,2)
            for iM = 1:3
                straind(iR,iC,iM,:) = [E11d(iR,iC,iM);   E22d(iR,iC,iM);    E33d(iR,iC,iM);...
                    2*E12d(iR,iC,iM); 2*E13d(iR,iC,iM);  2*E23d(iR,iC,iM)];
                if exist('xEBSD','var')
                    %solve the (HR-EBSD) boundary condition
                    K1=straind(iR,iC,iM,1)-straind(iR,iC,iM,3);
                    K2=straind(iR,iC,iM,2)-straind(iR,iC,iM,3);
                    K3=straind(iR,iC,iM,4)*Maps.Stiffness(4,3)+...
                        straind(iR,iC,iM,5)*Maps.Stiffness(5,3)+...
                        straind(iR,iC,iM,6)*Maps.Stiffness(6,3);
                    
                    e33n=-(K1*Maps.Stiffness(1,3)+K2*Maps.Stiffness(2,3)+K3)/...
                        (Maps.Stiffness(1,3)+Maps.Stiffness(2,3)+Maps.Stiffness(3,3));
                    e11n=K1+e33n;
                    e22n=K2+e33n;
                    %new strain/stress vector
                    co = straind(iR,iC,iM,4:6);
                    straind(iR,iC,iM,:)=[e11n;e22n;e33n;co(:)];
                end
                stressd(iR,iC,iM,:) = Maps.Stiffness*squeeze(straind(iR,iC,iM,:));
                sXXd(iR,iC,iM) = stressd(iR,iC,iM,1);
                sYYd(iR,iC,iM) = stressd(iR,iC,iM,2);
                sZZd(iR,iC,iM) = stressd(iR,iC,iM,3);
                sXYd(iR,iC,iM) = stressd(iR,iC,iM,4);
                sXZd(iR,iC,iM) = stressd(iR,iC,iM,5);
                sYZd(iR,iC,iM) = stressd(iR,iC,iM,6);
            end
        end
    end
    E11d = squeeze(straind(:,:,:,1));
    E22d = squeeze(straind(:,:,:,2));
    E33d = squeeze(straind(:,:,:,3));
    E12d = 0.5*squeeze(straind(:,:,:,4));
    E13d = 0.5*squeeze(straind(:,:,:,5));
    E23d = 0.5*squeeze(straind(:,:,:,6));
else    % Chauchy stress tensor, assming linear-elastic, isotropic material for
    % plane stress conditions
    Maps.G = Maps.E/(2*(1 + Maps.nu));
    sXXd = Maps.E/(1-Maps.nu^2)*(E11d+Maps.nu*(E22d+E33d));
    sYYd = Maps.E/(1-Maps.nu^2)*(E22d+Maps.nu*(E11d+E33d));
    sZZd = Maps.E/(1-Maps.nu^2)*(E33d+Maps.nu*(E11d+E22d));
    sXYd = 2*Maps.G*E12d;
    sXZd = 2*Maps.G*E13d;
    sYZd = 2*Maps.G*E23d;
end
Wd = 0.5*(E11d.*sXXd+E22d.*sYYd+E33d.*sZZd+2*E12d.*sXYd+2*E13d.*sXZd+2*E23d.*sYZd);

plotDecomposedStrain(E11d,E22d,E33d,E12d,E13d,E23d,Maps);
if isfield(Maps,'SavingD')
saveas(gcf, [fileparts(Maps.SavingD) '\Strain_I_II_III.fig']);
saveas(gcf, [fileparts(Maps.SavingD) '\Strain_I_II_III.tif']);  close
end

%%
% Generate q field
celw = 1;  % Width of area contour. Has to be an odd number.
dQdX = ones(DataSize)/(Maps.stepsize*celw);
dQdX = flipud(tril(flipud(tril(dQdX))))-flipud(triu(flipud(triu(dQdX))));
dQdY = dQdX';
dQdX(dQdX.*dQdY~=0) = 0;
% Domain integral
dA = ones(DataSize).*Maps.stepsize^2;
JAd = ((sXXd.*uXXd+sXYd.*uYXd+sXZd.*uZXd-Wd).*dQdX+(sYYd.*uYXd+sXYd.*uXXd+...
    sYZd.*uZXd).*dQdY).*dA;
% Contour selection
mid = floor(DataSize(1)/2);
[a,b] = meshgrid(1:DataSize(1));
linecon=round(max((abs(a-mid-1/2)),abs(b-mid-1/2)));
% Summation of integrals
for ii = 1:mid-1
    areaID = linecon>=(ii-floor(celw/2)) & linecon<=(ii+floor(celw/2));
    J.JKIII(:,ii,:) = sum(sum(JAd.*areaID,'omitnan'),'omitnan');
end

%% Equivalent SIF
if exist('xEBSD','var')% Mode I & II plane strain conditions
    %     if strcmp(Maps.stressstate,'plane_stress')
    %         J.KRaw(1:2,:) = sqrt(abs(J.JKIII(1:2,:))*Maps.E);
    %     else
    J.KRaw(1:2,:) = sqrt(abs(J.JKIII(1:2,:))*Maps.E)/sqrt((1-Maps.nu^2));
    %     end
else
    J.KRaw(1:2,:) = sqrt(abs(J.JKIII(1:2,:))*Maps.E);
end

J.KRaw(3,:) = sqrt(J.JKIII(3,:)*2*Maps.G);      % Mode III
J.K.Raw = sum(abs(J.JKIII));

%%
figure; plot(J.Raw); hold; plot(J.K.Raw); legend('J','J_K')%trim acess
set(gcf,'position',[98 311 1481 667])
text(1:length(J.Raw),J.Raw,string([1:length(J.Raw)]))
oh = input('where to cut the contour? '); close

%%
J.Raw = J.Raw(1:oh);
J.K.Raw = J.K.Raw(1:oh);
Kd = sqrt(real(J.KRaw).^2+imag(J.KRaw).^2);
Kd = Kd(:,1:oh);
KI.Raw   = Kd(1,:)*1e-6;
KII.Raw  = Kd(2,:)*1e-6;
KIII.Raw = Kd(3,:)*1e-6;
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<1;       dic = 1;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div  = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);
KI.true  = round(mean(rmoutliers(real(KI.Raw(contrs:end)))),dic);
KI.div   = round(std(rmoutliers(real(KI.Raw(contrs:end))),1),dic);
KII.true = round(mean(rmoutliers(real(KII.Raw(contrs:end)))),dic);
KII.div  = round(std(rmoutliers(real(KII.Raw(contrs:end))),1),dic);
KIII.true= round(mean(rmoutliers(real(KIII.Raw(contrs:end)))),dic);
KIII.div = round(std(rmoutliers(real(KIII.Raw(contrs:end))),1),dic);

plotJKIII(KI,KII,KIII,J,Maps.stepsize/saf,Maps.units.xy)
if isfield(Maps,'SavingD')
    saveas(gcf, [fileparts(Maps.SavingD) '\J_KI_II_III_v2.fig']);
    saveas(gcf, [fileparts(Maps.SavingD) '\J_KI_II_III_v2.tif']);  %  close all
    save([fileparts(Maps.SavingD) '\KIII_2D_v2.mat'],'Maps','J','KI','KII','KIII','saf');
end
end

function [ alldata,dataum ] = reshapeStrain( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
x  = raw_data(:,1);
y  = raw_data(:,2);
z  = raw_data(:,3);
xVec = unique(x);
yVec = unique(y);
zVec = unique(z);
E11 = raw_data(:,4);
E12 = raw_data(:,5);
E13 = raw_data(:,6);
E21 = raw_data(:,7);
E22 = raw_data(:,8);
E23 = raw_data(:,9);
E31 = raw_data(:,10);
E32 = raw_data(:,11);
E33 = raw_data(:,12);

[dataum.X,dataum.Y,dataum.Z] = meshgrid(xVec,yVec,zVec);
[nRows, nCols , nDep] = size(dataum.X);
dataum.E11 = zeros(nRows, nCols, nDep); %Initialise
dataum.E12 = zeros(nRows, nCols, nDep); %Initialise
dataum.E13 = zeros(nRows, nCols, nDep); %Initialise
dataum.E21 = zeros(nRows, nCols, nDep); %Initialise
dataum.E22 = zeros(nRows, nCols, nDep); %Initialise
dataum.E23 = zeros(nRows, nCols, nDep); %Initialise
dataum.E31 = zeros(nRows, nCols, nDep); %Initialise
dataum.E32 = zeros(nRows, nCols, nDep); %Initialise
dataum.E33 = zeros(nRows, nCols, nDep); %Initialise

for iRow = 1:nRows % loop rows
    for iCol = 1:nCols % loop cols
        for iDep = 1:nDep
            xt = dataum.X(iRow,iCol,iDep);
            yt = dataum.Y(iRow,iCol,iDep);
            zt = dataum.Z(iRow,iCol,iDep);
            idx = find(x==xt & y==yt & z==zt);
            if ~isempty(idx)
                E11t = E11(idx(1));
                E12t = E12(idx(1));
                E13t = E13(idx(1));
                E21t = E21(idx(1));
                E22t = E22(idx(1));
                E23t = E23(idx(1));
                E31t = E31(idx(1));
                E32t = E32(idx(1));
                E33t = E33(idx(1));
                dataum.E11(iRow,iCol,iDep) = E11t;
                dataum.E12(iRow,iCol,iDep) = E12t;
                dataum.E13(iRow,iCol,iDep) = E13t;
                dataum.E21(iRow,iCol,iDep) = E21t;
                dataum.E22(iRow,iCol,iDep) = E22t;
                dataum.E23(iRow,iCol,iDep) = E23t;
                dataum.E31(iRow,iCol,iDep) = E31t;
                dataum.E32(iRow,iCol,iDep) = E32t;
                dataum.E23(iRow,iCol,iDep) = E33t;
            end
        end
    end
end
alldata = [dataum.X(:)      dataum.Y(:)     dataum.Z(:)     dataum.E11(:) ...
           dataum.E12(:)    dataum.E13(:)  	dataum.E21(:)   dataum.E22(:)   ...
           dataum.E23(:)    dataum.E31(:)   dataum.E32(:)   dataum.E33(:)]; 

end

function [E,v,G,Co] = effectiveE_nu(C) 
% this function caclulates effective Youn modulus and Possion ratio for a
% ansitropic material based on this paper
% Reference: https://doi.org/10.3390/cryst8080307
BV = (C(1,1)+2*C(1,2))/3;               % Voigt bulk modulus
GV = (C(1,1)-C(1,2)+3*C(4,4))/5;        % Voigt shear modulus

S = C^-1;
BR = 1/(3*S(1,1)+6*S(1,2));             % Reuss bulk modulus
GR = 5/(4*S(1,1)-4*S(1,2)+3*S(4,4));   % Reuss shear modulus

B = (BR+BV)/2;                          % Hill’s average bulk modulus
G = (GR+GV)/2;                          % Hill’s average shear modulus
E = 9*B*G/(3*B+G);                      % Young’s modulus (E)
v = (3*B-E)/(6*B);                      % Poisson’s ratio
Co = [];

% K = (C(1,1)+C(2,2)+C(3,3)+2*(C(1,2)+C(2,3)+C(1,2)))/9; % istropic shear Modulus
% Gv = (C(1,1)+C(2,2)+C(3,3)-(C(1,2)+C(2,3)+C(1,2))+2*(C(4,4)+C(5,5)+C(6,6)))/15; % Bulk Modulus

%% Paper: What is the Young’s Modulus of Silicon?
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
if G<0 || E<0 || v<0 
    E = (Co.Ex+Co.Ey)/2;
    v = (Co.vxz+Co.vyz)/2;
    G = (Co.Gxz+Co.Gyz)/2;
end
end

function [Crop] = CroppingEqually(Maps)
close all;                  fig=subplot(1,1,1);
imagesc(Maps.X(1,:),Maps.Y(:,1),Maps.E11);
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
    eval(sprintf('Crop.E%d%d = Maps.E%d%d(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));',...
        iV,iO,iV,iO));
    end
end

%% XY, steps and stifness
Crop.X   = Maps.X(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.Y   = Maps.Y(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.Z   = Maps.Z(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.X   = Crop.X - min(min(Crop.X));  	Crop.Y   = Crop.Y - min(min(Crop.Y));
if (Crop.X(1) - Crop.X(end))>0;         Crop.X   = flip(Crop.X,2);         end
if (Crop.Y(1) - Crop.Y(end))>0;         Crop.Y   = flip(Crop.Y,1);         end

end

