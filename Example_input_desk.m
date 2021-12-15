% Example, synthetic data
[MatProp,~,alldata] = Calibration_2DKIII(3,1,2);

[J,KI,KII,KIII] = KIII_2D(alldata,MatProp); % as desigignated maps
[J,KI,KII,KIII] = KIII_2D(MatProp);         % or just MatProp

%% HR-EBSD data
filename = 'S3_4_20kV_15nA_200ms_XEBSD';
[Maps,alldata]=GetGrainData(filename);
KIII_2D(alldata,Maps);