%%  Test for separability,
%       using the 3*3 bound entangled states by P. Horodecki
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: Gilbert_SLOCC.m
%   Copyright @ 2017 Jiangwei Shang and Otfried Guehne
%   Email: jiangwei.shang@quantumlah.org
%   Last updated: August 06, 2017

%%
% ****************************SLOCC STRUCTURES*****************************
%
dimHL = 3; % local Hilbert space dimension
nPart = 2; % number of parties
npure = 1; % number of nonequivalent (under SLOCC) pure states 

% global Hilbert space dimension
dges = dimHL^nPart;
% The shuffle mapping used to change to different parties in the SLOCC optimization
iMap = ShuffleInit(nPart,dimHL);
% Allocate the PsiDecomp array containing the pure states
PsiDecomp(1:dges,1:npure) = zeros;


% ******************DEFINE DECOMPOSITION STATES****************************
% The set of pure states in the decomposition which are nonequivalent under SLOCC


% %%(i) Product state for any multipartite system: Fully separable testing
% %%|00..0>
PsiDecomp(1:dges,1) = zeros;
PsiDecomp(1,1) = 1.;


% % %%(ii) For testing biseparability we define the three biseparable bell states
% %%|ABC>
% phiAB(1:8) = zeros;
% phiAB(1) = 1./sqrt(2.);
% phiAB(7) = 1./sqrt(2.);
% phiAC(1:8) = zeros;
% phiAC(1) = 1./sqrt(2.);
% phiAC(6) = 1./sqrt(2.);
% phiBC(1:8) = zeros;
% phiBC(1) = 1./sqrt(2.);
% phiBC(4) = 1./sqrt(2.);
% 
% PsiDecomp(1:dges,1:1) = phiAB;
% PsiDecomp(1:dges,2:2) = phiAC;
% PsiDecomp(1:dges,3:3) = phiBC;


% % %%(iii) Testing biseparability of 4 qubit states needs more states
% % 0
% prod(1:2) = zeros;
% prod(1) = 1;
% 
% % phi+ Bell state
% Bell(1:4) = zeros;
% Bell(1) = 1/sqrt(2);
% Bell(4) = 1/sqrt(2);
% 
% GHZ(1:8) = zeros;
% % 000
% GHZ(1) = 1/sqrt(2);
% % 111
% GHZ(8) = 1/sqrt(2);
% 
% dum(1:dges,1:dges) = zeros;
% % A-BCD
% PsiDecomp(1:dges,1:1) = (kron(prod,GHZ));
% % B-CDA
% [PsiDecomp(1:dges,2:2),dum] = ShuffleR(PsiDecomp(1:dges,1:1),dum,iMap,dges);
% % C-DAB
% [PsiDecomp(1:dges,3:3),dum] = ShuffleR(PsiDecomp(1:dges,2:2),dum,iMap,dges);
% % D-ABC
% [PsiDecomp(1:dges,4:4),dum] = ShuffleR(PsiDecomp(1:dges,3:3),dum,iMap,dges);
% % AB-CD
% PsiDecomp(1:dges,5:5) = (kron(Bell,Bell));
% % BC-DA
% [PsiDecomp(1:dges,6:6),dum]=ShuffleR(PsiDecomp(1:dges,5:5),dum,iMap,dges);
% % AC-BD 
% PsiDecomp(1:dges,7:7) = zeros;
% PsiDecomp(1,7:7) = 1/2;
% PsiDecomp(6,7:7) = 1/2;
% PsiDecomp(11,7:7) = 1/2;
% PsiDecomp(16,7:7) = 1/2;


%%
% for demonstration
dim = 9;
id = eye(dim)/dim;
r_sep = 1/8; % radius of the sep ball in dimension 3*3

% The family of 3*3 bound entangled states by P. Horodecki
a = 0.3;
rho_PH = 1/(8*a+1)*[a 0 0 0 a 0 0 0 a;
                 0 a 0 0 0 0 0 0 0;
                 0 0 a 0 0 0 0 0 0;
                 0 0 0 a 0 0 0 0 0;
                 a 0 0 0 a 0 0 0 a;
                 0 0 0 0 0 a 0 0 0;
                 0 0 0 0 0 0 (1+a)/2 0 sqrt(1-a^2)/2;
                 0 0 0 0 0 0 0 a 0;
                 a 0 0 0 a 0 sqrt(1-a^2)/2 0 (1+a)/2];

             
p = 0.86;             
rho = p*rho_PH + (1-p)*id;
l2 = 0.5*sqrt(sum(svd(rho-id).^2)); % Hilbert-Schmidt norm

% SLOCC structures under consideration
opts = struct;
opts.dimHL = dimHL;
opts.nPart = nPart;
opts.npure = npure;
opts.PsiDecomp = PsiDecomp; 
% Basic settings for Gilbert_SLOCC.m
opts.imaxG = 1000; % maximum number of iterations
opts.tol = 1e-5; % tolerance: distance threshold
opts.mm = 50; % the memory

epsilon = 1e-3;
for i = 1:100
    rho_t = (1+epsilon)*rho - epsilon*id;
    l1 = 0.5*sqrt(sum(svd(rho-rho_t).^2)); % Hilbert-Schmidt norm
    [bitt, ~, dist] = Gilbert_SLOCC(rho_t, opts);
    
    if bitt==1
        r_t = dist*l2/l1;
        if r_t <= r_sep
            fprintf('rho is separable.\n');
            break;
        elseif epsilon > 2*opts.tol
            epsilon = epsilon/2;
        else
            fprintf('Cannot decide!\n');
            break;
        end
    else
        fprintf('Gilbert algo did not converge, consider increasing the number of iter.\n');
    end

end

