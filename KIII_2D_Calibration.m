function [alldata,MatProp] = KIII_2D_Calibration(KI,KII,KIII)
sz = 100;

stressstat = 'plane_stress';
MatProp.E = 210e9;                                % Young's Modulus
MatProp.nu = 0.3;                                  % Poisson ratio
    switch stressstat
        case 'plane_strain' % for xEBSD
            kappa = 3 - (4 .* MatProp.nu); % [/]
        case 'plane_stress'
            kappa = (3 - MatProp.nu)./(1 + MatProp.nu); % [/]
    end                                                           % Bulk modulus
G = MatProp.E/(2*(1 + MatProp.nu));                                                          % Shear modulus
% SIF loading
KI = KI*1e6;                                                                     % Mode I SIF
KII = KII*1e6;                                                                    % Mode II SIF
KIII = KIII*1e6;     
% Maps = [1/Maps.E          -Maps.nu/Maps.E     -Maps.nu/Maps.E 0 0 0
%                  -Maps.nu/Maps.E        1/Maps.E        -Maps.nu/Maps.E 0 0 0
%                  -Maps.nu/Maps.E    -Maps.nu/Maps.E         1/Maps.E    0 0 0
%                   0 0 0     2*(1+Maps.nu)/Maps.E                          0 0
%                   0 0 0 0       2*(1+Maps.nu)/Maps.E                        0
%                   0 0 0 0 0         2*(1+Maps.nu)/Maps.E];
% Maps = Maps^-1;

%% Anayltical displacement data.
stepsize = 1/sz*2;
lin = stepsize*(ceil(-1/stepsize)+1/2):stepsize:stepsize*(floor(1/stepsize)-1/2);
[X,Y,Z] = meshgrid(lin,lin,0);
[th,r] = cart2pol(X,Y);

% displacement data
Ux = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))));
Uy = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))));
Uz = ( 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2));

%% JMAN approach (without FEM) - Standard J-integral.
[E11,E12,E13] = crackgradient(Ux,stepsize);
[E21,E22,E23] = crackgradient(Uy,stepsize);
[E31,E32,E33] = crackgradient(Uz,stepsize);
alldata = [X(:) Y(:) Z(:) E11(:) E12(:) E13(:) E21(:) E22(:) E23(:) E31(:) E32(:) E33(:)]; 
end
%% Support function
function [cx,cy,cz]=crackgradient(c,dx)
c=squeeze(c);
[row,~]=size(c);
midr=floor(row/2);
ctop=c(1:midr,:);
cbot=c(midr+1:end,:);
[cxtop,cytop]=gradient(ctop,dx);
[cxbot,cybot]=gradient(cbot,dx);
cx=[cxtop;cxbot];
cy=[cytop;cybot];
cz=zeros(size(cx));
end