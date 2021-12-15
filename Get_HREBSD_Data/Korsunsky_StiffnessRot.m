function DA=Korsunsky_StiffnessRot(R,L)
%% Salvati, E., Sui, T. and Korsunsky, A.M., 2016.
% Uncertainty quantification of residual stress evaluation by the FIB-DIC
% ring-core method due to elastic anisotropy effects.
% International Journal of Solids and Structures, 87, pp.61-69. 
R11 = R(1,1);   R12 = R(1,2);   R13 = R(1,3);
R21 = R(2,1);   R22 = R(2,2);   R23 = R(2,3);
R31 = R(3,1);   R32 = R(3,2);   R33 = R(3,3);
TsA = [   R11^2,   R21^2,   R31^2,         2*R21*R31,         2*R11*R31,         2*R11*R21 
              R12^2,   R22^2,   R32^2,         2*R22*R32,         2*R12*R32,         2*R12*R22 
              R13^2,   R23^2,   R33^2,         2*R23*R33,         2*R13*R33,         2*R13*R23 
            R12*R13, R22*R23, R32*R33, R22*R33 + R23*R32, R12*R33 + R13*R32, R12*R23 + R13*R22 
            R11*R13, R21*R23, R31*R33, R21*R33 + R23*R31, R11*R33 + R13*R31, R11*R23 + R13*R21 
            R11*R12, R21*R22, R31*R32, R21*R32 + R22*R31, R11*R32 + R12*R31, R11*R22 + R12*R21];
        
TeA = [R11^2,     R21^2,     R31^2,           R21*R31,           R11*R31,           R11*R21
            R12^2,     R22^2,     R32^2,           R22*R32,           R12*R32,           R12*R22
            R13^2,     R23^2,     R33^2,           R23*R33,           R13*R33,           R13*R23
            2*R12*R13, 2*R22*R23, 2*R32*R33, R22*R33 + R23*R32, R12*R33 + R13*R32, R12*R23 + R13*R22
            2*R11*R13, 2*R21*R23, 2*R31*R33, R21*R33 + R23*R31, R11*R33 + R13*R31, R11*R23 + R13*R21
            2*R11*R12, 2*R21*R22, 2*R31*R32, R21*R32 + R22*R31, R11*R32 + R12*R31, R11*R22 + R12*R21];
 DA=inv(TsA)*L*TeA;
end