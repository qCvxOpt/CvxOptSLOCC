%%  dg_Gilbert
%      The dg_Gilbert algorithm
%
%   rhoS_opt: The optimal state returned
%   stats: Various statistics returned
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: Gilbert_SLOCC.m, qmt.m
%   Copyright @ 2017 Jiangwei Shang and Otfried Guehne
%   Email: jiangwei.shang@quantumlah.org
%   Last updated: August 06, 2017

function [rhoS_opt, stats] = dg_Gilbert(operators, f, opts)
%
% grab dimensions
if ~isa(operators,'cell') % nonseparable
    operators = {operators};
end

dims = zeros(numel(operators),2);
for i=1:numel(operators)
    temp = size(operators{i});
    dims(i,1) = temp(1);
    dims(i,2) = temp(end);
end
d = prod(dims(:,1));
K = prod(dims(:,2));
assert(K==numel(f));


% process options
defaults = struct;

% common options
% %------------ for Gilbert ------------%
defaults.rho_in = eye(d)/d;

% SLOCC structures under consideration
defaults.dimHL = 3; % local Hilbert space dimension
defaults.nPart = 2; % number of parties
defaults.npure = 1; % number of nonequivalent (under SLOCC) pure states
defaults.PsiDecomp = [1; zeros(d-1,1)]; % fully separable testing
% Basic settings for Gilbert_SLOCC.m
defaults.imaxG = 1000; % maximum number of iterations
defaults.tol = 1e-5; % tolerance: distance threshold
defaults.mm = 50; % the memory
% %------------ end Gilbert ------------%


defaults.imax = 100; % maximum number of iterations
defaults.statptprec = 1e-10; % precision
defaults.tol_eps = 1e-3; % tolerance for stepsize: epsilon

defaults.epsilon = 25; % initial epsilon
defaults.eps_factor = 0.9;


if ~exist('opts','var')
    opts = defaults;
else
    % scan for invalid options
    names = fieldnames(opts);
    for i=1:numel(names)
        if ~isfield(defaults, names{i})
            error('dg_Gilbert:opts','unknown option %s',names{i});
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


% %------------ for Gilbert_SLOCC.m ------------%
optsG = struct;
% SLOCC structures under consideration
optsG.dimHL = opts.dimHL; % local Hilbert space dimension
optsG.nPart = opts.nPart; % number of parties
optsG.npure = opts.npure; % number of nonequivalent (under SLOCC) pure states
optsG.PsiDecomp = opts.PsiDecomp; % fully separable testing
% Basic settings for Gilbert_SLOCC.m
optsG.imaxG = opts.imaxG; % maximum number of iterations
optsG.tol = opts.tol; % tolerance: distance threshold
optsG.mm = opts.mm; % the memory
% %-------------------- end --------------------%


% the output
stats = struct;
stats.fvalS = zeros(opts.imax,1); % for projected separable states
stats.fval = zeros(opts.imax,1); % for states to be projected
stats.dist = zeros(opts.imax,1);
stats.bitt = zeros(opts.imax,1); % converged or not for Gilbert
stats.eps_step = zeros(opts.imax,1); % record stepsize in each update


% the data
% f0 = f; % keep the original f
fmap = (f~=0);
f = f(fmap);

% initialization
rhoS = opts.rho_in;
rhoS_opt = rhoS; % the target

pbS = qmt(rhoS, operators);
stats.fvalS0 = -f'*log(pbS(fmap));
stats.best_fvalS = stats.fvalS0;

%
for i = 1:opts.imax
    
    % gradient on rhoS
    adj = zeros(K,1);
    adj(fmap) = f./pbS(fmap);
    rmatrix = qmt(adj, operators, 'adjoint'); 
    
    % to update
    while true
        % update rho, based on rhoS
        rho_temp = (eye(d)+opts.epsilon/2*(rmatrix-eye(d)))*rhoS*...
            (eye(d)+opts.epsilon/2*(rmatrix-eye(d)));
        rho = rho_temp/real(trace(rho_temp));

        pb = qmt(rho, operators);
        stats.fval(i) = -f'*log(pb(fmap)); % new fval value
        
        if stats.fval(i) > stats.best_fvalS
            fprintf('Walks too much!!!\n');
            opts.epsilon = opts.epsilon*opts.eps_factor; % reduce the stepsize
        else
            break; % to proceed
        end
    end
    stats.eps_step(i) = opts.epsilon;

    % to find the closest separable state rhoS w.r.t. rho
    [opts.bitt(i), rhoS, opts.dist(i)] = Gilbert_SLOCC(rho, optsG);

    pbS = qmt(rhoS, operators);
    stats.fvalS(i) = -f'*log(pbS(fmap)); % new fvalS value

    % to check
    diff_fvalS = stats.fvalS(i) - stats.best_fvalS;
    if diff_fvalS > 0 && opts.epsilon > opts.tol_eps
        opts.epsilon = opts.epsilon*opts.eps_factor; % reduce the stepsize
        rhoS = rhoS_opt;
        stats.fvalS(i) = stats.best_fvalS;
    elseif diff_fvalS <= 0 && opts.epsilon > opts.tol_eps % updated
        rhoS_opt = rhoS;
        stats.best_fvalS = stats.fvalS(i);
    else
        fprintf('Cannot update anymore!!!\n');
        break;
    end
    
    fprintf('iter = %d, epsilon = %f, fvalS = %f, diff_fvalS = %f\n\n',i,stats.eps_step(i),stats.fvalS(i),diff_fvalS);
    
    % threshold
    if opts.epsilon < opts.tol_eps %|| abs(diff_fvalS) < opts.statptprec
        break;
    end

end

stats.fvalS = stats.fvalS(1:i);
stats.fval = stats.fval(1:i);
stats.dist = stats.dist(1:i);
stats.bitt = stats.dist(1:i);
stats.eps_step = stats.eps_step(1:i);

end

