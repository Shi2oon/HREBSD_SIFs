function [Maps] = rotateStrains(Maps,L,B)
% this code take either theta (L) degree of rotation or the burger vector
% (B)
% and line vecttor (L) of slip band to reslove the stress compoenent on it
% the variable in should have all the components seperated

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

A = {Maps.A11 Maps.A12 Maps.A13;Maps.A21 Maps.A22 Maps.A23;Maps.A31 Maps.A23 Maps.A33};
[txf.A]    = componentTensorTransform(A,R);
Maps.A11 =  txf.A{1,1};    Maps.A12 =  txf.A{1,2};   	Maps.A13 =  txf.A{1,3};
Maps.A21 =  txf.A{2,1};  	Maps.A22 =  txf.A{2,2};    Maps.A23 =  txf.A{2,3};
Maps.A31 =  txf.A{3,1};    Maps.A32 =  txf.A{3,2};  	Maps.A33 =  txf.A{3,3};
%{
S = {Maps.S11 Maps.S12 Maps.S13;Maps.S21 Maps.S22 Maps.S23;Maps.S31 Maps.S23 Maps.S33};
[txf.S]    = componentTensorTransform(S,R);
Maps.S11 =  txf.S{1,1};	Maps.S12 =  txf.S{1,2};	Maps.S13 =  txf.S{1,3};
Maps.S21 =  txf.S{2,1};	Maps.S22 =  txf.S{2,2};	Maps.S23 =  txf.S{2,3};
Maps.S31 =  txf.S{3,1};	Maps.S32 =  txf.S{3,2};	Maps.S33 =  txf.S{3,3};

E = {Maps.E11 Maps.E12 Maps.E13;Maps.E21 Maps.E22 Maps.E23;Maps.E31 Maps.E23 Maps.E33};
[txf.E]    = componentTensorTransform(E,R);
Maps.E11 =  txf.E{1,1};   	Maps.E12 =  txf.E{1,2};    Maps.E13 =  txf.E{1,3}; 
Maps.E21 =  txf.E{2,1};   	Maps.E22 =  txf.E{2,2};    Maps.E23 =  txf.E{2,3};
Maps.E31 =  txf.E{3,1};   	Maps.E32 =  txf.E{3,2};    Maps.E33 =  txf.E{3,3};

W = {Maps.W11 Maps.W12 Maps.W13;Maps.W21 Maps.W22 Maps.W23;Maps.W31 Maps.W23 Maps.W33};
[txf.W]    = componentTensorTransform(W,R);
Maps.W11 =  txf.W{1,1};   	Maps.W12 =  txf.W{1,2};    Maps.W13 =  txf.W{1,3};       
Maps.W21 =  txf.W{2,1};   	Maps.W22 =  txf.W{2,2};    Maps.W23 =  txf.W{2,3}; 
Maps.W31 =  txf.W{3,1};   	Maps.W32 =  txf.W{3,2};    Maps.W33 =  txf.W{3,3}; 
%}
end

%%
