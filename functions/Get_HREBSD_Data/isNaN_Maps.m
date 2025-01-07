function Maps = isNaN_Maps(Maps)
answers = input('Do you want to trim out Grain boundries (Y/N)? ','s');
if answers == 'Y' || answers == 'y'
    DoS = MoreOutNaN(Maps.S12,unique(Maps.X));

Maps.GND(isnan(DoS))=NaN;      Maps.PH(isnan(DoS))=NaN;       Maps.MAE(isnan(DoS))=NaN;
% Rotation
Maps.W11(isnan(DoS))=NaN;      Maps.W12(isnan(DoS))=NaN;      Maps.W13(isnan(DoS))=NaN; 
Maps.W21(isnan(DoS))=NaN;      Maps.W22(isnan(DoS))=NaN;      Maps.W23(isnan(DoS))=NaN;  
Maps.W31(isnan(DoS))=NaN;      Maps.W32(isnan(DoS))=NaN;      Maps.W33(isnan(DoS))=NaN;  
% Stress
Maps.S11(isnan(DoS))=NaN;      Maps.S12(isnan(DoS))=NaN;      Maps.S13(isnan(DoS))=NaN; 
Maps.S21(isnan(DoS))=NaN;      Maps.S22(isnan(DoS))=NaN;      Maps.S23(isnan(DoS))=NaN;  
Maps.S31(isnan(DoS))=NaN;      Maps.S32(isnan(DoS))=NaN;      Maps.S33(isnan(DoS))=NaN; 
% Strain
Maps.E11(isnan(DoS))=NaN;      Maps.E12(isnan(DoS))=NaN;      Maps.E13(isnan(DoS))=NaN; 
Maps.E21(isnan(DoS))=NaN;      Maps.E22(isnan(DoS))=NaN;      Maps.E23(isnan(DoS))=NaN;  
Maps.E31(isnan(DoS))=NaN;      Maps.E32(isnan(DoS))=NaN;      Maps.E33(isnan(DoS))=NaN;
% Deformation gradient
Maps.A11(isnan(DoS))=NaN;      Maps.A12(isnan(DoS))=NaN;      Maps.A13(isnan(DoS))=NaN; 
Maps.A21(isnan(DoS))=NaN;      Maps.A22(isnan(DoS))=NaN;      Maps.A23(isnan(DoS))=NaN; 
Maps.A31(isnan(DoS))=NaN;      Maps.A32(isnan(DoS))=NaN;      Maps.A33(isnan(DoS))=NaN;
end