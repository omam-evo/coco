clear

% path to the MATLAB build
addpath('..\opt','..\opt\subroutines');

% only suit name (bbob), no logger
rng(1);
suite_name = 'bbob';
suite = cocoSuite(suite_name,'', '');
fcts = 24; inst = 15; dims = 6;

% select problem id, show what problem it is
get_f = 17;
get_N = 4;  %[2,3,5,10,20,40]
get_inst = 1;
BUDGET_MULTIPLIER = 1e4;

VERBOSE = 1;
STOP.G_MAX = 1e9;
STOP.X_TOL = 1e-12;
STOP.SIGMA_D_STOP = 1e-12;
STOP.COND_SQRT = 1e7;
STOP.SIGMA_STOP = nan;
STOP.F_STAG_TOL = 1e-12;

NUM_PLOTS_SHOW = 1;

problem_id = (get_inst-1)+(get_f-1)*inst + fcts*inst*(get_N-1); 
bbob = cocoSuiteGetProblem(suite, problem_id);
str_problem = cocoProblemGetName(bbob);
F_STOP = cocoProblemGetBestValue(bbob);
N = cocoProblemGetDimension(bbob);

%get dimension, lower bounds, upper bounds and a starting point
lb = cocoProblemGetSmallestValuesOfInterest(bbob);
ub = cocoProblemGetLargestValuesOfInterest(bbob);
lb = lb'; ub = ub';
x0 = cocoProblemGetInitialSolution(bbob); 

%define the function to minimize from the problem
fun = @(x) cocoEvaluateFunction(bbob,x);

%optimize with your favourite method
% options = optimset('Display','iter','MaxFunEvals',200*dim);
% [x_res,f_res] = fminsearch(fun,x0,options) %Nelder-Mead
%[x_res,f_res] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options) %BFGS
%[x_res,f_res] = myoptimizer(...) %favourite optimizer

i = 1;
% figure; hold on;
while BUDGET_MULTIPLIER*N > cocoProblemGetEvaluations(bbob)
    
    feval_remain = BUDGET_MULTIPLIER*N-cocoProblemGetEvaluations(bbob);
    STOP.FEVAL_MAX = feval_remain;

    %% Init
    lambda= 4+floor(3*log(double(N)));
    mu = floor(lambda/2); 
    x0 = rand(N,1).*8-4;
    sigma0 = 2;

    if i<=NUM_PLOTS_SHOW
        SHOW_PLOW = 1;
    else
        SHOW_PLOW = 0;
    end

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
    [x,f,~] =  apop_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, 1, VERBOSE);

    %% LRA-CMA-ES
    % [x,f,~] =  lra_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% ETA-CMA-ES
    % STOP.SIGMA_STOP = 1e-8;
    % lambda = max(10,double(N));
    % mu = lambda/2;
    % [x,f,~] =  eta_cma_es(fun, bbob, mu, lambda, x0, sigma0, STOP, VERBOSE);

    %% Apply offset of -F_STOP to obtain strictly positive f-dynamics 
    % => needed for logarithmic axis
    if i<=NUM_PLOTS_SHOW
        fig = gcf;
        lineObjs = findobj(fig, 'type', 'line');
        if strcmpi(lineObjs(end).DisplayName, 'f')
            lineObjs(end).YData = lineObjs(end).YData-F_STOP;
        end
    end

    if cocoProblemFinalTargetHit(bbob)
        break;
    end
    % plot(f_g-F_STOP, '-'); 
    % plot(sigD_g, 'k:'); 

    i = i+1;
end
%yscale('log');
%xlabel('g'); ylabel('f(g)-F');


if cocoProblemFinalTargetHit(bbob)
    disp(['success: ',str_problem]);
else
    disp(['fail: ',str_problem]);
end

%free the suit (otherwise, MATLAB crashes might happen)
cocoSuiteFree(suite);

%surf plot
% if dim == 2
%     x1_discretization = linspace(lb(1),ub(1));
%     x2_discretization = linspace(lb(2),ub(2));
%     [X1,X2] = meshgrid(x1_discretization,x2_discretization);
%     Z = 0*X1;
%     for it=1:numel(X1)
%         Z(it) = fun([X1(it), X2(it)]);
%     end
%     surf(X1,X2,Z);
% end
