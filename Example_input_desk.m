% Example, synthetic data
clc;clear;close all
[~,~,alldata] = Calibration_2DKIII(3,1,2);
               
Prop.E = 210e9; 
Prop.nu = 0.3;              
Prop.units.St = 'Pa';
Prop.units.xy = 'm';
Prop.stressstat = 'plane_stress';

[J,KI,KII,KIII] = KIII_2D_v2(alldata,Prop); % as desigignated maps
% [J,KI,KII,KIII] = KIII_2D(MatProp);         % or just MatProp

%% HR-EBSD data, works only with data extract from xEBSD Matlab code and you need to download MTEX up to version 5.8
clc;clear;close all
% restoredefaultpath;clc;clear;close all
% run('A:\OneDrive - Nexus365\GitHub\mtex-5.2.beta2\install_mtex.m')
filename = 'Crack_in_Si_XEBSD';
[Maps,alldata] = GetGrainData(filename);
KIII_2D_v2(alldata,Maps);
