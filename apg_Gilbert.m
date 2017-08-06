%%  apg_Gilbert
%      The apg_Gilbert algorithm
%
%   rho: The optimal state returned
%   stats: Various statistics returned
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: Gilbert_SLOCC.m, qmt.m
%   Copyright @ 2017 Jiangwei Shang and Otfried Guehne
%   Email: jiangwei.shang@quantumlah.org
%   Last updated: August 06, 2017

function [rho, stats] = apg_Gilbert(operators, f, opts)
%
% APG_GILBERT Optimization over the SLOCC region using APG
%
%   RHO = QT_APG(OPERATORS, F)
%
%   RHO = QT_APG(OPERATORS, F, OPTS)
%
%   [RHO, STATS] = QT_APG(...)

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
defaults.rho0 = 'Gilbert'; % initial density matrix
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


defaults.threshold_step = 0.5*sqrt(d)*d*eps; % stop when trace distance upper bound between current and previous iterate is less than this
defaults.threshold_fval = -f(f~=0)'*log(f(f~=0)); % stop when negative log likeilhood less than this
defaults.threshold_dist = 0; % stop when trace distance between rho_star and rho is less than or equal to this
defaults.rho_star = [];
defaults.imax = 100; % maximum number of iterations

% instrumentation
defaults.save_rhos = false;
defaults.save_cond_stats = false;

% tweak behavior type of algorithm
defaults.guard = true; % guard against varrho exiting
defaults.fix_gradient = false; % prevent giant components of the gradient
defaults.bb = true; % use Barzilai-Bowein to initialize step sizes instead of afactor increases (seems to work well! at least more stable)
defaults.late_restart = false; % set to true to keep last step when restarting
defaults.keep_theta = false; % retain value of theta after a restart
defaults.traceless = false; % make the gradient traceless
defaults.precond = false; % enable experimental preconditioner
defaults.restart = 'fval'; % grad vs. fval restarting
defaults.accel = 'fista_tfocs'; % 'fista_tfocs','fista','apg' different ways to grow and compute beta

% tweak behavior parameters 
defaults.t0 = []; % set to a scalar if we know a good starting t, set to [] for Lipschitz estimation
defaults.restart_grad_param = 0.01; % how much extra sensitivity (>=0.0, <1.0)
defaults.bfactor = 0.5; % t backtracks by this much
defaults.afactor = 1.1; % t grows by this much
defaults.minimum_t = eps; % stop if t gets smaller than this
defaults.bootstrap_threshold = 0.01; % threshold for bootstrap using CG (smaller value -> longer to run CG)

if ~exist('opts','var')
    opts = defaults;
else
    % scan for invalid options
    names = fieldnames(opts);
    for i=1:numel(names)
        if ~isfield(defaults, names{i})
            error('apg_Gilbert:opts','unknown option %s',names{i});
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


% discard zero-valued frequencies (they have no effect on the merit function)
% f0 = f; % keep the original f
fmap = (f~=0);
f = f(fmap);
coeff_temp = zeros(K,1);

start = tic;

stats = struct;
stats.steps = zeros(opts.imax,1); % upper bound on trace distance between current iterate and previous iterate
stats.fvals = zeros(opts.imax,1);
stats.dists = zeros(opts.imax,1); % lower bound or actual trace distance between current iterate and rho_star
stats.dists_true = false(opts.imax,1); % set to true if actual trace distance
stats.comp_prob = 0;
stats.comp_grad = 0;
stats.comp_fval = 0;
stats.comp_proj = 0;
stats.times = zeros(opts.imax,1);
stats.thetas = zeros(opts.imax,1);
stats.ts = zeros(opts.imax,1);
stats.satisfied_step = false;
stats.satisfied_fval = false;
stats.satisfied_dist = false;


if opts.save_rhos
    stats.rhos = {};
end

if ischar(opts.rho0)
    switch(opts.rho0)
        case 'Gilbert'
            rho = opts.rho_in;
            % Initial varrho is the MLE, not used
            % Initial rho is varrho_projected
        case 'white'
            rho = eye(d)/d;
        case 'bootstrap'
            % run CG with line search until "condition number" of
            % adjustment vector stops changing by much
            opts2=struct; opts2.imax = opts.imax; opts2.mincondchange=opts.bootstrap_threshold; 
            opts2.threshold_fval = opts.threshold_fval;
            opts2.threshold_dist = opts.threshold_dist;
            opts2.rho_star = opts.rho_star;
            opts2.save_rhos = opts.save_rhos;
            coeff_temp(fmap) = f;
            time_offset = toc(start);
            [rho, boot_stats] = CGLS(operators,coeff_temp,opts2);
            boot_stats.times = boot_stats.times + time_offset;
            opts.imax = opts.imax-numel(boot_stats.fvals); % reduce number of iterations we can run
            if boot_stats.satisfied_fval || boot_stats.satisfied_dist
                % CG run was good enough
                stats = boot_stats;
                return;
            end
        otherwise
            error('qt_apg:initializer',['unknown initializer specification: ',opts.rho0]);
    end
else
    rho = opts.rho0;
end

% set varrho to rho in the initial step
varrho = rho; 

varrho_changed = true;
gradient = [];
theta = 1;
t = opts.t0;
bb_okay = false; % don't do two-step Barzilai-Borwein on first step
precond = opts.precond;

% compute initial probabilities
probs_rho = qmt(rho, operators);
probs_rho = probs_rho(fmap);

% record the fval from rho_ML projected
stats.fval0 = -f'*log(probs_rho);
fval = stats.fval0;

stats.comp_prob = stats.comp_prob+1;
stats.comp_fval = stats.comp_fval+1;
probs_varrho = probs_rho;

stats.best_rho = rho;
stats.best_fval = fval;

if opts.save_cond_stats
    hessian_proxy = [];
end

for i=1:opts.imax
    % compute new gradient if varrho has changed
    if varrho_changed
        if opts.bb && i>1
            old_gradient = gradient;
        end
        coeff_temp(fmap) = f./probs_varrho;
        if opts.save_cond_stats 
            if ~isempty(hessian_proxy)
                old_hessian_proxy = hessian_proxy;
                hessian_proxy = f./(probs_varrho).^2;
                stats.cond_change_angle(i,1) = real(acos(real(old_hessian_proxy(:)'*hessian_proxy(:))/norm(old_hessian_proxy(:))/norm(hessian_proxy(:))));
            else
                hessian_proxy = f./(probs_varrho).^2;
            end
        end

        if opts.fix_gradient
            % fix any ridiculously strong gradient components
            kappa = (max(f)/min(f));
            xi = kappa./(-f.*log(f/kappa));
            coeff_temp(fmap) = coeff_temp(fmap).*(-xi.*f.*log(probs_varrho))./(coeff_temp(fmap) - xi.*f.*log(probs_varrho));
        end
        
        gradient = -qmt(coeff_temp, operators, 'adjoint');
        if opts.traceless
            gradient = gradient - trace(gradient)*eye(d)/d;
            gradient = gradient - trace(gradient)*eye(d)/d;
        end
        stats.comp_grad = stats.comp_grad+1;
        fval_varrho = -f'*log(probs_varrho);
        stats.comp_fval = stats.comp_fval+1;
        if precond
            g = gradient;
            gradient = gradient*varrho' + varrho*gradient';
        end
        if opts.bb && bb_okay
            varrho_diff = varrho(:)-old_varrho(:);
            gradient_diff = gradient(:)-old_gradient(:);
            denominator = gradient_diff'*gradient_diff;
            if denominator > 0
                t = abs(real(varrho_diff'*gradient_diff))/denominator;
            end
        end
    end
    
    if i==1 && exist('boot_stats','var')
        boot_stats.fvals(end) = fval_varrho;
    end
    
    if isempty(t)
        % compute local Lipschitz constant from derivatives
        probs_gradient = qmt(-gradient, operators);
        probs_gradient = probs_gradient(fmap);
        first_deriv = -f'*(probs_gradient./probs_varrho);
        second_deriv = f'*(probs_gradient.^2./(probs_varrho.^2));
        stats.comp_prob = stats.comp_prob+1;
        t = -first_deriv/second_deriv;
    else
        if ~(opts.bb && bb_okay)
            t = t * opts.afactor;
        end
    end
    
    % backtrack for finding a good value of t
    t_good = false;
    fval_new = [];
    t_threshold = opts.minimum_t / norm(gradient(:));
    while ~t_good
        if ~isempty(fval_new)
            if precond
                gradient = g;
                precond = false;
            else
            new_t_estimate = second_order/max(0,fval_new-fval_varrho-first_order);
            % this comparison is false if any term is NaN (namely
            % new_t_estimate)
            if new_t_estimate > t_threshold
                t = min(t*opts.bfactor, new_t_estimate);
            else
                % if new_t_estimate is NaN or less than or equal to t_threshold
                t = t*opts.bfactor;
            end
            end
        end        
        
        [~, rho_new, ~] = Gilbert_SLOCC(varrho - t*gradient, optsG);
        
        stats.comp_proj = stats.comp_proj + 1;
        probs_rho_new = qmt(rho_new, operators);
        probs_rho_new = probs_rho_new(fmap); 
        fval_new = -f'*log(probs_rho_new);
        stats.comp_prob = stats.comp_prob+1;
        stats.comp_fval = stats.comp_fval+1;
        stats.fvals(i) = fval_new;
        delta = rho_new(:) - varrho(:);
        first_order = real(gradient(:)'*delta(:));
        second_order = 0.5*delta(:)'*delta(:);
        % multiplied by t so that we don't get NaN or Inf if t is too small
        t_good = ~(t*fval_new > t*fval_varrho + t*first_order + second_order); % not greater than catches NaNs
    end
    
    if fval_new < stats.best_fval
        stats.best_fval = fval_new;
        stats.best_rho = rho_new;
    end

    if opts.save_rhos
        stats.rhos{i,1} = rho_new;
    end
    
    stats.ts(i) = t;
    stats.thetas(i) = theta;
    
    fprintf('iter = %d, fvalS = %f\n\n',i,stats.best_fval);

    % check threshold
    stats.steps(i) = 0.5*sqrt(d)*norm(rho_new-rho,'fro');
    stats.satisfied_step = stats.steps(i) <= opts.threshold_step;
    stats.satisfied_fval = fval_new <= opts.threshold_fval;
    if ~isempty(opts.rho_star)
        stats.dists(i) = 0.5*norm(rho_new-opts.rho_star,'fro');
        if stats.dists(i) <= opts.threshold_dist
            % do additional check with actual trace distance
            stats.dists(i) = 0.5*sum(svd(rho_new-opts.rho_star));
            stats.dists_true(i) = true;
            stats.satisfied_dist = stats.dists(i) <= opts.threshold_dist;
        end
    end
    if t < t_threshold || stats.satisfied_step || stats.satisfied_fval || stats.satisfied_dist
        stats.times(i) = toc(start);
        rho = rho_new;
        break;
    end

    % record previous value of varrho for Barzilai-Borwein
    if opts.bb
        old_varrho = varrho;
    end
    
    % adaptive restart stuff
    switch(opts.restart)
        case 'fval'
            do_restart = (fval_new > fval);
        case 'grad'
            vec1 = varrho(:)-rho_new(:);
            vec2 = rho_new(:)-rho(:);
            do_restart = real(vec1'*vec2) > -opts.restart_grad_param*norm(vec1)*norm(vec2);
        case 'none'
            do_restart = false;
        otherwise
            error('qt_apg:restarttype','erroneous restart type');
    end

    bb_okay = ~do_restart; % enable Barzilai-Borwein if no restart  
    
    if do_restart
        if opts.late_restart
            % keep the last step
            rho = rho_new;
            probs_rho = probs_rho_new;
            fval = fval_new;
        end
        varrho = rho;
        probs_varrho = probs_rho;
        varrho_changed = (theta>1);
        if ~opts.keep_theta
            theta = 1;
        end
        stats.times(i) = toc(start);
        stats.fvals(i) = fval;
        if opts.save_rhos
            stats.rhos{i,1} = rho;
        end
        continue;
    end
        
    % acceleration
    if i>1 && stats.ts(i) > eps
        Lfactor = stats.ts(i-1)/stats.ts(i);
    else
        Lfactor = 1;
    end
    switch(opts.accel)
        case 'fista'
            theta_new = (1+sqrt(1+4*theta^2))/2;
            beta = (theta-1)/theta_new;
            theta = theta_new;
        case 'fista_tfocs'
            theta_hat = sqrt(Lfactor)*theta;
            theta_new = (1+sqrt(1+4*theta_hat^2))/2;
            beta = (theta_hat-1)/theta_new;
            theta = theta_new;
        case 'apg' 
            theta_new = theta*(sqrt(1+theta^2/4)-theta/2);
            beta = theta*(1-theta)/(theta^2+theta_new);
            theta = theta_new;
        case 'none'
            beta = 0;
        otherwise
            error('qt_apg:accel','unknown acceleration scheme');
    end
    
    % update    
    varrho = rho_new + beta*(rho_new-rho);
    probs_varrho = probs_rho_new + beta*(probs_rho_new-probs_rho);
    if opts.guard && min(probs_varrho) <= 0
        % discard momentum if momentum causes varrho to become infeasible
        % retain theta to keep estimate of current condition number
        varrho = rho_new;
        probs_varrho = probs_rho_new;
    end
    varrho_changed = true;
    rho = rho_new;
    probs_rho = probs_rho_new;
    fval = fval_new;
    stats.times(i) = toc(start);
end

stats.fvals = stats.fvals(1:i);
stats.steps = stats.steps(1:i);
stats.dists = stats.dists(1:i);
stats.dists_true = stats.dists_true(1:i);
stats.times = stats.times(1:i);
stats.thetas = stats.thetas(1:i);
stats.ts = stats.ts(1:i);

if exist('boot_stats','var')
    if opts.save_rhos
        stats.rhos = [boot_stats.rhos;stats.rhos];
    end
    if opts.save_cond_stats
        stats.cond_change_angle = [boot_stats.cond_change_angle;stats.cond_change_angle];
    end
    stats.fvals = [boot_stats.fvals;stats.fvals];
    stats.steps = [boot_stats.steps;stats.steps];
    stats.dists = [boot_stats.dists;stats.dists];
    stats.dists_true = [boot_stats.dists_true;stats.dists_true];
    stats.times = [boot_stats.times;stats.times];
    stats.comp_prob = stats.comp_prob + boot_stats.comp_prob;
    stats.comp_grad = stats.comp_grad + boot_stats.comp_grad;
    stats.comp_fval = stats.comp_fval + boot_stats.comp_fval;
    stats.comp_proj = stats.comp_proj + boot_stats.comp_proj;
end

end
