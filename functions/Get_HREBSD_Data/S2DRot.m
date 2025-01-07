function S_trans = S2DRot(L,theta)
%% Salvati, E., Sui, T. and Korsunsky, A.M., 2016.
% Uncertainty quantification of residual stress evaluation by the FIB-DIC
% ring-core method due to elastic anisotropy effects.
% International Journal of Solids and Structures, 87, pp.61-69. 
R(1,1) = cosd(theta);
R(1,2) = sind(theta);
R(2,1) = -sind(theta);
R(2,2) = cosd(theta);
R(3,3) = 1;
S_trans=Korsunsky_StiffnessRot(R,L);
end
% end
% plot(strike, c22, '.-')
% title('Anisotropy Stiffness Element C22')
% xlabel('\theta^o')
% ylabel('Stiffness (GPa)')
% xlim([0, max(strike)])