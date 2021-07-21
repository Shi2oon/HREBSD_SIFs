function [Maps] = rotateStrains(Maps,L,B)
% this code take either theta (L) degree of rotation or the burger vector
% (B)
% and line vecttor (L) of slip band to reslove the stress compoenent on it
% the variable in should have all the components seperated
S = {Maps.S11 Maps.S12 Maps.S13;Maps.S21 Maps.S22 Maps.S23;Maps.S31 Maps.S23 Maps.S33};
E = {Maps.E11 Maps.E12 Maps.E13;Maps.E21 Maps.E22 Maps.E23;Maps.E31 Maps.E23 Maps.E33};
W = {Maps.W11 Maps.W12 Maps.W13;Maps.W21 Maps.W22 Maps.W23;Maps.W31 Maps.W23 Maps.W33};
A = {Maps.A11 Maps.A12 Maps.A13;Maps.A21 Maps.A22 Maps.A23;Maps.A31 Maps.A23 Maps.A33};
%{
Strain(:,:,1,1) = Maps.E11;	Strain(:,:,1,2) = Maps.E12;	Strain(:,:,1,3) = Maps.E13;
Strain(:,:,2,1) = Maps.E21;	Strain(:,:,2,2) = Maps.E22;	Strain(:,:,2,3) = Maps.E23;
Strain(:,:,3,1) = Maps.E31;	Strain(:,:,3,2) = Maps.E32;	Strain(:,:,3,3) = Maps.E33;

Stress(:,:,1,1) = Maps.S11;	Stress(:,:,1,2) = Maps.S12;	Stress(:,:,1,3) = Maps.S13;
Stress(:,:,2,1) = Maps.S21;	Stress(:,:,2,2) = Maps.S22;	Stress(:,:,2,3) = Maps.S23;
Stress(:,:,3,1) = Maps.S31;	Stress(:,:,3,2) = Maps.S32;	Stress(:,:,3,3) = Maps.S33;

RotW(:,:,1,1) = Maps.W11;	RotW(:,:,1,2) = Maps.W12;	RotW(:,:,1,3) = Maps.W13;
RotW(:,:,2,1) = Maps.W21;	RotW(:,:,2,2) = Maps.W22;	RotW(:,:,2,3) = Maps.W23;
RotW(:,:,3,1) = Maps.W31;	RotW(:,:,3,2) = Maps.W32;	RotW(:,:,3,3) = Maps.W33;
%}
if exist('B','var')
    N = cross(B,L);
    R = [  B(1) L(1) N(1)
        B(2) L(2) N(2)
        B(3) L(3) N(3)];
elseif size(L,2)==3
    R = L; % rotation matrix
else
    Rx = [1 0 0; 0 cosd(L) -sind(L); 0 sind(L) cosd(L)];
    %     Ry = [cosd(L) 0 sind(L); 0 1 0; -sind(L) 0 cosd(L)];
    %     Rz = [cosd(L) -sind(L) 0; sind(L) cosd(L) 0; 0 0 1];
    R = transpose(Rx);
end
[txf.A]    = componentTensorTransform(A,R);
[txf.E]    = componentTensorTransform(E,R);
[txf.S]    = componentTensorTransform(S,R);
[txf.W]    = componentTensorTransform(W,R);
Maps.A11 =  txf.A{1,1};    Maps.A12 =  txf.A{1,2};   	Maps.A13 =  txf.A{1,3};
Maps.A21 =  txf.A{2,1};  	Maps.A22 =  txf.A{2,2};    Maps.A23 =  txf.A{2,3};
Maps.A31 =  txf.A{3,1};    Maps.A32 =  txf.A{3,2};  	Maps.A33 =  txf.A{3,3};

Maps.S11 =  txf.S{1,1};	Maps.S12 =  txf.S{1,2};	Maps.S13 =  txf.S{1,3};
Maps.S21 =  txf.S{2,1};	Maps.S22 =  txf.S{2,2};	Maps.S23 =  txf.S{2,3};
Maps.S31 =  txf.S{3,1};	Maps.S32 =  txf.S{3,2};	Maps.S33 =  txf.S{3,3};

Maps.E11 =  txf.E{1,1};   	Maps.E12 =  txf.E{1,2};    Maps.E13 =  txf.E{1,3}; 
Maps.E21 =  txf.E{2,1};   	Maps.E22 =  txf.E{2,2};    Maps.E23 =  txf.E{2,3};
Maps.E31 =  txf.E{3,1};   	Maps.E32 =  txf.E{3,2};    Maps.E33 =  txf.E{3,3};

Maps.W11 =  txf.W{1,1};   	Maps.W12 =  txf.W{1,2};    Maps.W13 =  txf.W{1,3};       
Maps.W21 =  txf.W{2,1};   	Maps.W22 =  txf.W{2,2};    Maps.W23 =  txf.W{2,3}; 
Maps.W31 =  txf.W{3,1};   	Maps.W32 =  txf.W{3,2};    Maps.W33 =  txf.W{3,3}; 
    %{
    for iC=1:size(Strain,1)
        for iR=1:size(Strain,2)
            newStrain(iC,iR,:,:) = R*squeeze(Strain(iC,iR,:,:))*transpose(R);
            newStress(iC,iR,:,:) = R*squeeze(Stress(iC,iR,:,:))*transpose(R);
            newRotW(iC,iR,:,:)   = R*squeeze(RotW(iC,iR,:,:))*transpose(R);
        end
    end
    
Maps.E11 = newStrain(:,:,1,1); Maps.E12 = newStrain(:,:,1,2); Maps.E13 = newStrain(:,:,1,3);
Maps.E21 = newStrain(:,:,2,1); Maps.E22 = newStrain(:,:,2,2); Maps.E23 = newStrain(:,:,2,3);
Maps.E31 = newStrain(:,:,3,1); Maps.E32 = newStrain(:,:,3,2); Maps.E33 = newStrain(:,:,3,3);

Maps.S11 = newStress(:,:,1,1); Maps.S12 = newStress(:,:,1,2); Maps.S13 = newStress(:,:,1,3);
Maps.S21 = newStress(:,:,2,1); Maps.S22 = newStress(:,:,2,2); Maps.S23 = newStress(:,:,2,3);
Maps.S31 = newStress(:,:,3,1); Maps.S32 = newStress(:,:,3,2); Maps.S33 = newStress(:,:,3,3);

Maps.W11 = newRotW(:,:,1,1); Maps.W12 = newRotW(:,:,1,2); Maps.W13 = newRotW(:,:,1,3);
Maps.W21 = newRotW(:,:,2,1); Maps.W22 = newRotW(:,:,2,2); Maps.W23 = newRotW(:,:,2,3);
Maps.W31 = newRotW(:,:,3,1); Maps.W32 = newRotW(:,:,3,2); Maps.W33 = newRotW(:,:,3,3);
    %}
end