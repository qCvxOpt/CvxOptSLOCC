%%  Gilbert_SLOCC
%      The Gilbert_SLOCC algorithm
%
%   bitt = 1: converged; bitt = 0: not converged
%   s_cls: closest SLOCC state w.r.t. r
%   dist_cls: HS distance btwn s_cls and r
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: ShuffleInit.m, RandomLocalUnitary.m, WSLOCCMaxG.m
%   Copyright @ 2017 Jiangwei Shang and Otfried Guehne
%   Email: jiangwei.shang@quantumlah.org
%   Last updated: August 06, 2017

function [bitt, s_cls, dist_cls] = Gilbert_SLOCC(r, opts)
%
dim = size(r,1); % global dimension

% process options
defaults = struct;

% SLOCC structures under consideration
defaults.dimHL = 3; % local Hilbert space dimension
defaults.nPart = 2; % number of parties
defaults.npure = 1; % number of nonequivalent (under SLOCC) pure states
defaults.PsiDecomp = [1; zeros(dim-1,1)]; % fully separable testing
% Basic settings for Gilbert_SLOCC.m
defaults.imaxG = 1000; % maximum number of iterations
defaults.tol = 1e-5; % tolerance: distance threshold
defaults.mm = 50; % the memory

if ~exist('opts','var')
    opts = defaults;
else
    % scan for invalid options
    names = fieldnames(opts);
    for i=1:numel(names)
        if ~isfield(defaults, names{i})
            error('Gilbert_SLOCC:opts','unknown option %s',names{i});
        end
    end
    % populate defaults if not in opts
    names = fieldnames(defaults);
    for i=1:numel(names)
        if ~isfield(opts, names{i})
            opts.(names{i}) = defaults.(names{i});
        end
    end
end

% the shuffle mapping used to change to different parties in the SLOCC optimization
iMap = ShuffleInit(opts.nPart,opts.dimHL);

Lopt = zeros(dim,dim,opts.mm);
A = zeros(dim*dim,opts.mm+1);
dist = zeros(1,opts.imaxG);

%
% -------------------------------- %
% building up the initial memory
s = eye(dim)/dim; % step 1
for i = 1:opts.mm
    
    % step 2: the maximal overlap with a product vector
    iil = mod(i,opts.npure)+1;
    WPrimeD = RandomLocalUnitary(opts.PsiDecomp(1:dim,iil),opts.nPart,opts.dimHL);
    WPrime = WSLOCCMaxG(r-s+2*eye(dim),WPrimeD,iMap,opts.dimHL,opts.nPart); 
    pv = WPrime*WPrime';

    Lopt(:,:,i) = pv; % memory added in

    % step 3: to update
    temp_eps = real(trace((r-s)*(pv-s)))/real(trace((pv-s)*(pv-s)));
    eps_k = min(temp_eps,1);
    s = (1-eps_k)*s + eps_k*pv;
    
    % fprintf('mm %d: dist = %f\n',i,0.5*sum(svd(r-s)));
end

%
% the main loop
bitt = 0;
dist(1) = 0.5*sqrt(sum(svd(r-s).^2)); % Hilbert-Schmidt norm
dist_cls = dist(1);
s_cls = s;
for i = 1:opts.imaxG
    
    fprintf('Step %d: dist = %f; diff = %f\n',i,dist(i),dist_cls-dist(i));
    
    % stopping criterion
    if dist_cls < opts.tol
        bitt = 1;
        break;
    end
      
    % step 2: the maximal overlap with a product vector
    iil = mod(i,opts.npure)+1;
    WPrimeD = RandomLocalUnitary(opts.PsiDecomp(1:dim,iil),opts.nPart,opts.dimHL);
    WPrime = WSLOCCMaxG(r-s+2*eye(dim),WPrimeD,iMap,opts.dimHL,opts.nPart); 
    pv = WPrime*WPrime';

    % the memory
    for j = 1:opts.mm-1
        Lopt(:,:,j) = Lopt(:,:,j+1);
        A(:,j) = reshape(Lopt(:,:,j),[dim*dim, 1]);
    end
    Lopt(:,:,opts.mm) = pv;
    A(:,opts.mm) = pv(:);
    A(:,opts.mm+1) = s(:);

    % solve linear least squares
    lb = zeros(opts.mm+1,1);
    ub = ones(opts.mm+1,1);
    Aeq = ones(1,opts.mm+1);
    beq = 1;
    options = optimset('LargeScale','off','Display','off');
    x = real(lsqlin(A,r(:),[],[],Aeq,beq,lb,ub,[],options));

    % step 3: to update
    s = reshape(A*x, [dim, dim]);
    dist(i+1) = 0.5*sqrt(sum(svd(r-s).^2)); % Hilbert-Schmidt norm
    
    if dist(i+1) <= dist_cls
        dist_cls = dist(i+1);
        s_cls = s;
    end
    
end

end
