%%  Demonstration of the likelihood-ratio test (LRT), 
%       using the 3*3 bound entangled states by P. Horodecki
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: Gilbert_SLOCC.m
%   Copyright @ 2017 Jiangwei Shang and Otfried Guehne
%   Email: jiangwei.shang@quantumlah.org
%   Last updated: August 06, 2017

%%
% SIC-POVM for qutrit
d3 = 3;
w3 = exp(1j*2*pi/d3);
t3 = exp(1j*pi*(d3+1)/d3);

X3 = [0 0 1; 1 0 0; 0 1 0];
Z3 = [1 0 0; 0 w3 0; 0 0 w3'];

% fid in d=3
fid3 = 1/sqrt(2)*[0; 1; -1];

% to build
D3 = zeros(d3,d3,d3^2);
SIC3_kets = zeros(d3,d3^2);
SIC3 = zeros(d3,d3,d3^2);
for i = 1:d3^2
    ind1 = ceil(i/d3);
    ind2 = 1+rem(i-1,d3);
    D3(:,:,i) = Z3^ind1*X3^ind2*t3^(ind1*ind2);
    SIC3_kets(:,i) = D3(:,:,i)*fid3/sqrt(d3);
    SIC3(:,:,i) = SIC3_kets(:,i)*SIC3_kets(:,i)';
end

POVM = repmat({SIC3_kets},[1 2]); % for 2-qutrit

%%
% basic settings
dim = 9;
id = eye(dim)/dim;

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
             
p = 1.;
rho_t = p*rho_PH + (1-p)*id; % the true state
prob_t = qmt(rho_t,POVM);

%
% Quantum State Tomography
rho_ML = rho_t; % set the true state as MLE

% the data D
f = prob_t; % set f==p
fmap = (f~=0);

probs_ML = qmt(rho_ML, POVM);
fval_ML = -f(fmap)'*log(probs_ML(fmap));

%%
% %---------------------------- dg_Gilbert ----------------------------% %
%
% SLOCC structures under consideration
opts1 = struct;
opts1.dimHL = 3; % local Hilbert space dimension
opts1.nPart = 2; % number of parties
opts1.npure = 1; % number of nonequivalent (under SLOCC) pure states
opts1.PsiDecomp = [1; zeros(dim-1,1)]; % fully separable testing
% Basic settings for Gilbert_SLOCC.m
opts1.imaxG = 100; % maximum number of iterations
opts1.tol = 1e-5; % tolerance: distance threshold
opts1.mm = 50; % the memory

opts1.rho_in = id;
opts1.imax = 50;
[rho_dg, stats_dg] = dg_Gilbert(POVM, f, opts1);


%%
% %---------------------------- apg_Gilbert ----------------------------% %
%
% SLOCC structures under consideration
opts2 = struct;
opts2.dimHL = 3; % local Hilbert space dimension
opts2.nPart = 2; % number of parties
opts2.npure = 1; % number of nonequivalent (under SLOCC) pure states
opts2.PsiDecomp = [1; zeros(dim-1,1)]; % fully separable testing
% Basic settings for Gilbert_SLOCC.m
opts2.imaxG = 20; % maximum number of iterations
opts2.tol = 1e-5; % tolerance: distance threshold
opts2.mm = 50; % the memory

opts2.rho_in = id;
opts2.imax = 50;
opts2.accel = 'fista_tfocs'; % 'fista_tfocs','fista','apg' for beta
[rho_apg, stats_apg] = apg_Gilbert(POVM, f, opts2);


%%
% The likelihood-ratio test (LRT)
% lrt = 2*(F_S-F_ML), as lrt=-2*log(L_S/L_ML)

lrt_dg = 2*(stats_dg.fvalS-fval_ML);
lrt_apg = 2*(stats_apg.fvals-fval_ML);

% to plot
semilogy(lrt_dg,'.-')
hold on
semilogy(lrt_apg,'x-r')

set(gca,'fontsize',18,'FontName','Times New Roman')

%axis([0 60 0.0306 0.0338]);
xlabel('Steps','FontSize',22,'FontName','Times New Roman');
ylabel('$\lambda/N$','Interpreter','latex','FontSize',22,'FontName','Times New Roman');
legend('DG','APG','Location','NorthWest','FontSize',20,'FontName','Times New Roman');

