clear
tic

%% Generic
rng(1);
addpath('subroutines')

%% Config
VERBOSE = 1;
STOP.FEVAL_MAX = 1e7;
STOP.G_MAX = 1e9;       % use for noisy sphere
STOP.X_TOL = 1e-12;
STOP.SIGMA_D_STOP = 1e-12;
STOP.COND_SQRT = 1e7;
STOP.SIGMA_STOP = nan;
STOP.F_STAG_TOL = 1e-12;

%% Problem
TRIALS = 5;
% FIT = Fitness('sphere', []); N = 40;
% FIT = Fitness('Ellipsoid-H', []); N = 40;

FIT = Fitness('rastrigin', [10,2*pi]); N = 20;
% FIT = Fitness('rastrigin', [1,2*pi]); N = 100;

% FIT = Fitness('noisysphere', []); STOP.G_MAX = 1e4;  N=40;

fun = @(x) FIT.get_f(x);
STOP.F_STOP = 1e-8;

%% BBOB benchmarking?
bbob = [];
% lb = -5*ones(1,N);
% ub = 5*ones(1,N);

%% ES parameters
lambda = 4+floor(3*log(N)); %4+floor(3*log(N));
mu = floor(lambda/2);

n_success = 0;
for i = 1:TRIALS

    fprintf('main_dyn: %i/%i \n',i,TRIALS);
    % y = rand(N,1).*(ub-lb)' + lb';  
    x0 = 100*ones(N,1); %10*rand(N,1);
    lb = -5*ones(N,1);
    ub = 5*ones(N,1);
    sigma0 = 100; %rand*9+1;  % U~[1,10]
 
    %% (mu/mu_I,lambda)-sigmaSA-ES
    % [x,f,~] = sa_es(fun, bbob, mu, lambda, x0, sigma0, 1/sqrt(2*N), STOP, VERBOSE);

    %% purecmaes basic from Hansen
    % [x,f,~] = purecmaes_basic(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% BIPOP-CMA-ES
    % [x,f,~] =  bipop_cma_es(fun, bbob, mu, lambda, x0, lb, ub, sigma0, STOP, VERBOSE);

    %% CMSA-ES
    % [x,f,~] =  cmsa_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% MA-ES
    % ma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% Meta-ES
    % isolation_factor = 10;
    % lambda = 128;
    % mu = lambda/2;
    % [x,f,~] = meta_es(fun,bbob, mu, lambda, isolation_factor, x0, sigma0, STOP, VERBOSE);

    %% PSA-CMA-ES
    % [x,f,~] =  psa_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% APOP-CMA-ES
    % [x,f,~] =  apop_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, 1, VERBOSE);

    %% LRA-CMA-ES
    % [x,f,~] =  lra_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% ETA-CMA-ES
    % STOP.SIGMA_STOP = 1e-8;
    lambda = max(10,N);
    mu = lambda/2;
    [x,f,~] =  eta_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    if f <= STOP.F_STOP % not needed for coco
        n_success = n_success + 1;
    end
end % repeat problem

% myfigsize(gcf,8,5,9,9,0.5);
P_S = n_success/TRIALS;