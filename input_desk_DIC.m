clc;clear;close all
% Make sure the grid spacing is equal in x and y direction, also make sure
% the crack is at the centre of the map
Prop.E = 210e9;  % Young's Modulus, you can also input a stifness matrix
Prop.nu = 0.3;     % Poisson ratio         
Prop.units.St = 'Pa';  
Prop.units.xy = 'mm';% meter (m) or milmeter (mm) or micrometer(um);
Prop.stressstat = 'plane_stress';% 'plane_stress' OR 'plane_strain'
Prop.Operation = 'DIC'; % 'DIC' for raw 2D and 3D DIC data
DataDirect = fullfile(pwd,'Data','1KI-2KII-3KII_Data.dat'); % file locatio
Prop.results = fullfile(pwd,'Data','1KI-2KII-3KII_Data');

% end of user input
Data = importdata(DataDirect);
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(Data.data,Prop);