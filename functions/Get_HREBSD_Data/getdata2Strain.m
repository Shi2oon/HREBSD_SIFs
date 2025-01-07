function [Maps,alldata] = getdata2Strain(alldata,Maps,SavingD)
%% Prep Data % Strain   
[ ~,dataum ] = reshapeData( alldata(:,1:4));
    Maps.E11 = dataum.Ux;            Maps.E22 = dataum.Uy;
    [ ~,dataum ] = reshapeData([alldata(:,1:3) alldata(:,5)]);
    Maps.E12 = dataum.Uy;       Maps.E21 = Maps.E12; 
    
    E = zeros(size(Maps.E11,1),size(Maps.E11,2),3,3);
    Maps.E13 = E(:,:,2,2);      Maps.E23 = E(:,:,2,3);
    Maps.E31 = E(:,:,3,1);      Maps.E32 = E(:,:,3,2);      Maps.E33 = E(:,:,3,3);
    
    % Rotation
    W = zeros(size(Maps.E11,1),size(Maps.E11,2),3,3);
    Maps.W11 = W(:,:,1,1);      Maps.W12 = W(:,:,1,2);      Maps.W13 = W(:,:,1,3);
    Maps.W21 = W(:,:,2,1);      Maps.W22 = W(:,:,2,2);      Maps.W23 = W(:,:,2,3);
    Maps.W31 = W(:,:,3,1);      Maps.W32 = W(:,:,3,2);      Maps.W33 = W(:,:,3,3);
    % Stress
    S = zeros(size(Maps.E11,1),size(Maps.E11,2),3,3);
    Maps.S11 = S(:,:,1,1);      Maps.S12 = S(:,:,1,2);      Maps.S13 = S(:,:,1,3);
    Maps.S21 = S(:,:,2,1);      Maps.S22 = S(:,:,2,2);      Maps.S23 = S(:,:,2,3);
    Maps.S31 = S(:,:,3,1);      Maps.S32 = S(:,:,3,2);      Maps.S33 = S(:,:,3,3);

    % Deformation gradient
    A = zeros(size(Maps.E11,1),size(Maps.E11,2),3,3);
    Maps.A11 = A(:,:,1,1);      Maps.A12 = A(:,:,1,2);      Maps.A13 = A(:,:,1,3);
    Maps.A21 = A(:,:,2,1);      Maps.A22 = A(:,:,2,2);      Maps.A23 = A(:,:,2,3);
    Maps.A31 = A(:,:,3,1);      Maps.A32 = A(:,:,3,2);      Maps.A33 = A(:,:,3,3);
    
    Maps.GND    = zeros(size(squeeze(S(:,:,1,1))));
    Maps.RefID = Maps.A31;   Maps.PH = Maps.A31;     Maps.MAE = Maps.A31;   
    
    % stifness:  crystal orientation is defined as the rotation that transforms crystal
    % coordinates, i.e., a description of a vector or a tensor with respect to the crystal
    % reference frame, into specimen coordinates, i.e., a desciption of the same object
    % with respect to a specimen fixed reference frame.
    if Maps.type == 'A'
       Maps.nu  =  Maps.Stiffness(1,2)/(Maps.Stiffness(1,1)+ Maps.Stiffness(1,2));
       Maps.E   =  Maps.Stiffness(1,1)*(1-2*Maps.nu)*(1+Maps.nu)/(1-Maps.nu);
    else
        Maps.Stiffness = Maps.E*ones(6,6);
    end
    % Dim
    Maps.X   = dataum.X1;    Maps.Y   = dataum.Y1;
    Maps.stepsize  =(abs(Maps.X(1,1)-Maps.X(1,2)));
    Maps.Wo  = (1/2).*(Maps.S11.*Maps.E11 + Maps.S12.*Maps.E12 + Maps.S13.*Maps.E13 +...
        Maps.S21.*Maps.E21 + Maps.S22.*Maps.E22 + Maps.S23.*Maps.E23 +...
        Maps.S31.*Maps.E31 + Maps.S32.*Maps.E32 + Maps.S33.*Maps.E33);
    % units (defualt xEBSD units)
    Maps.units.xy = Maps.input_unit;      	Maps.units.S  = 'GPa';
    Maps.units.E  = 'Abs.';                 Maps.units.W = 'rad';

%% crop data
answer = input('Do you want to crop (C) and/or rotate (R) data (C/R/N)? ','s');
if answer == 'C' || answer == 'c'
    [Maps] = Cropping(Maps,SavingD);
    SavingD = [SavingD '\Cropped Data.mat'];
elseif answer == 'R' || answer == 'r'
    [Maps] = Rot2Crop(Maps,SavingD,answer);
    SavingD = [SavingD '\Crop & Rot Data.mat'];
else
    SavingD = [SavingD '\Full Data.mat'];
end

%% crack cordinates
% try;    Maps = isNaN_Maps(Maps);       end % trim GB
close all;
imagesc(Maps.X(1,:),Maps.Y(:,1),Maps.E12);
c=colorbar; c.Label.String = ['E_{12} [' Maps.units.xy ']'];
%     caxis([-5e-3 5e-3]);
set(gca,'Ydir','normal');	axis image;colormap jet
title('Answer in the command line');
xlabel(['X [' Maps.units.xy ' ]'],'FontSize',20,'FontName','Times New Roman');
ylabel(['Y [' Maps.units.xy ' ]'],'FontSize',20,'FontName','Times New Roman');
set(gcf,'WindowStyle','normal')
set(gcf,'position',[30 50 1300 950]);
title('E_{12} :: Select the Crack, start from crack tip');
[Maps.xo,Maps.yo] = ginput(2);
title('E_{12} :: Select the Crack mask, start from crack tip');
[Maps.xm,Maps.ym] = ginput(2); close

%% Get stifness tensor
Maps.SavingD = SavingD;
if strcmpi(answer, 'R')
    Maps.Stiffness  = S2DRot(Maps.Stiffness,Maps.theta); % rotate the stifness
end

%% save an exit
alldata = [Maps.X(:)   Maps.Y(:)   Maps.E11(:)   Maps.E22(:) Maps.E12(:)];
save(SavingD,'Maps','alldata'); % save