%% HR-EBSD data
clc;clear;close all
% restoredefaultpath;clc;clear;close all
% run('A:\OneDrive - Nexus365\GitHub\mtex-5.2.beta2\install_mtex.m')
filename = 'Crack_in_Si_XEBSD';
[Maps,alldata] = GetGrainData(filename);
KIII_2D_v2(alldata,Maps);

