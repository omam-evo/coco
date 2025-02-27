%%%%%%%%%%%%%%%%%%%%%%1e-12
% Generic Parameters %
%%%%%%%%%%%%%%%%%%%%%%
clear
addpath('../opt/subroutines','../opt');      % ES algorithms in "opt" folder
rng(1);                 % default rng

%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%
SET_ALGO = ['bipop'];         % PURE,APOP,BIPOP,PSA,LRA,META,ETA
labelrun = [SET_ALGO,'_l1']; % append suffix to label experiment
BUDGET_MULTIPLIER = 1e5;     % algorithm runs for BUDGET_MULTIPLIER*dimension funevals

STOP.G_MAX = 1e9;          % max number of generations (usually very large because function eval. are relevant)
STOP.X_TOL = 1e-12;        % min change in search space
STOP.SIGMA_D_STOP = 1e-12; % min product of sqrt. of largest eigenvalue (cov) and sigma
STOP.COND_SQRT = 1e7;      % max sqrt of covariance conditioning
STOP.SIGMA_STOP = nan;     % min sigma (can be set, but usually covered by SIGMA_D_STOP)
STOP.F_STAG_TOL = 1e-12;   % stagnation threshold of fitness values

%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Experiment    %
%%%%%%%%%%%%%%%%%%%%%%%%%

% choose a test suite and a matching logger, for
% example one of the following:
%
% bbob               24 unconstrained noiseless single-objective functions
% bbob-biobj         55 unconstrained noiseless bi-objective functions
% [bbob-biobj-ext*   92 unconstrained noiseless bi-objective functions]
% bbob-largescale    24 unconstrained noiseless single-objective functions in large dimensions
% bbob-constrained   48 constrained noiseless single-objective functions
% bbob-mixint        24 unconstrained noiseless single-objective functions with mixed-integer variables
% bbob-biobj-mixint  92 unconstrained noiseless bi-objective functions with mixed-integer variables
% bbob-boxed*         24 bound-constrained noiseless single-objective functions
% => suites with a star are partly implemented but may not yet be fully supported.
suite_name = 'bbob-largescale';
observer_name = 'bbob';
observer_options = strcat('result_folder: ',labelrun,...
    [' algorithm_name: RS '...
    ' algorithm_info: A_simple_random_search ']);

% initialize suite and observer with default options,
% to change the default, see 
% http://numbbo.github.io/coco-doc/C/#suite-parameters and
% http://numbbo.github.io/coco-doc/C/#observer-parameters
% for details.
suite = cocoSuite(suite_name, '', '');
observer = cocoObserver(observer_name, observer_options);

% set log level depending on how much output you want to see, e.g. 'warning'
% for fewer output than 'info'.
cocoSetLogLevel('info');
copyfile([mfilename,'.m'],fullfile('exdata',labelrun));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Experiment        %
%%%%%%%%%%%%%%%%%%%%%%%%%
tic
while true
    % get next problem and dimension from the chosen suite:
    bbob_problem = cocoSuiteGetNextProblem(suite, observer);
    if ~cocoProblemIsValid(bbob_problem)
        break;
    end
    fun = @(x) cocoEvaluateFunction(bbob_problem,x);
    N = double(cocoProblemGetDimension(bbob_problem));
    lb = cocoProblemGetSmallestValuesOfInterest(bbob_problem)';
    ub = cocoProblemGetLargestValuesOfInterest(bbob_problem)';
    str_problem = cocoProblemGetName(bbob_problem);
    disp(str_problem);

    i = 0; % count number of independent restarts
    delta = ub - lb;

    while BUDGET_MULTIPLIER*N > cocoProblemGetEvaluations(bbob_problem) + cocoProblemGetEvaluationsConstraints(bbob_problem)
        
        cocoObserverSignalRestart(observer, bbob_problem) % signal that a restart took place => why needed?
        
        feval_before = cocoProblemGetEvaluations(bbob_problem) + cocoProblemGetEvaluationsConstraints(bbob_problem);
        feval_remain = BUDGET_MULTIPLIER*N - feval_before;
        STOP.FEVAL_MAX = feval_remain;

        %% Init for all
        x0 = rand(N,1).*8-4;
        sigma0 = 2;

        if strcmpi(SET_ALGO, 'PURE')
            lambda= 4+floor(3*log(N)); mu = floor(lambda/2); 
            purecmaes_basic(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, 0);

        elseif strcmpi(SET_ALGO, 'PSA')
            lambda= 4+floor(3*log(N)); mu = floor(lambda/2); 
            psa_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, 0);

        elseif strcmpi(SET_ALGO, 'LRA')
            lambda= 4+floor(3*log(N)); mu = floor(lambda/2); 
            lra_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, 0);

        elseif strcmpi(SET_ALGO, 'APOP')
            map = containers.Map({'2','3','5','10','20','40'},[10,20,30,40,50,60]);
            lambda = map(num2str(N));
            %lambda= 4+floor(3*log(N)); 
            mu = floor(lambda/2); 
            APOP_MODE = 1;
            apop_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, APOP_MODE, 0);

            % NOTE: APOP in [NH17] uses default CMA-ES for initial run to gain performance on unimodal functions
            % if i==0  % first run simple CMA-ES ("init")
            %    APOP_MODE = 0;
            %    apop_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, APOP_MODE, 0);
            % else % ("reg")
            %    APOP_MODE = 1;
            %    apop_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, APOP_MODE, 0);
            % end

        elseif strcmpi(SET_ALGO, 'BIPOP')
            lambda= 4+floor(3*log(N)); mu = floor(lambda/2); 
            bipop_cma_es(fun, bbob_problem, mu, lambda, x0, lb, ub, sigma0, STOP, 0);

        elseif strcmpi(SET_ALGO, 'META')
            lambda = 32;
            mu = lambda/2;
            meta_es(fun, bbob_problem, mu, lambda, 10, x0, sigma0, STOP, 0);

        elseif strcmpi(SET_ALGO, 'ETA')
            STOP.SIGMA_STOP = 1e-8;
            lambda = max(10,N);
            mu = lambda/2;
            eta_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, 0);
        else
            error('Algorithm not found.')
        end

        % check whether things went wrong or whether experiment is over:
        feval_after = cocoProblemGetEvaluations(bbob_problem) + cocoProblemGetEvaluationsConstraints(bbob_problem);
        if cocoProblemFinalTargetHit(bbob_problem) == 1 || feval_after >= BUDGET_MULTIPLIER * N
            break;
        end
        if (feval_after == feval_before)
            fprintf('WARNING: Budget has not been used (%d/%d evaluations done)!\n', feval_before, BUDGET_MULTIPLIER * N);
            break;
        end
        if (feval_after < feval_before)
            fprintf('ERROR: Something weird happened here which should not happen: f-evaluations decreased');
        end
        i = i+1;
    end

end
cocoObserverFree(observer);
cocoSuiteFree(suite);
toc
