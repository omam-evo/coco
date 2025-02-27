clear

% path to the MATLAB build
addpath('..\mat_v4','..\mat_v4\es_for_meta');

% only suit name (bbob), no logger
rng(1);
suite_name = 'bbob';
suite = cocoSuite(suite_name,'', '');
fcts = 24; inst = 15; dims = 6;

% select problem id, show what problem it is
get_f = 17;
get_N = 4;  %[2,3,5,10,20,40]
get_inst = 1;
BUDGET_MULTIPLIER = 1e5;

STOP.G_MAX = 1e9;
STOP.X_TOL = 1e-12;
STOP.SIGMA_D_STOP = 1e-12;
STOP.COND_SQRT = 1e7;
STOP.SIGMA_STOP = nan;
STOP.F_STAG_TOL = 1e-12;

NUM_PLOTS_SHOW = 0;

problem_id = (get_inst-1)+(get_f-1)*inst + fcts*inst*(get_N-1); 
bbob_problem = cocoSuiteGetProblem(suite, problem_id);
str_problem = cocoProblemGetName(bbob_problem);
F_STOP = cocoProblemGetBestValue(bbob_problem);
N = cocoProblemGetDimension(bbob_problem);

%get dimension, lower bounds, upper bounds and a starting point
lb = cocoProblemGetSmallestValuesOfInterest(bbob_problem);
ub = cocoProblemGetLargestValuesOfInterest(bbob_problem);
lb = lb'; ub = ub';
x0 = cocoProblemGetInitialSolution(bbob_problem); 

%define the function to minimize from the problem
fun = @(x) cocoEvaluateFunction(bbob_problem,x);

%optimize with your favourite method
% options = optimset('Display','iter','MaxFunEvals',200*dim);
% [x_res,f_res] = fminsearch(fun,x0,options) %Nelder-Mead
%[x_res,f_res] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options) %BFGS
%[x_res,f_res] = myoptimizer(...) %favourite optimizer

i = 1;
% figure; hold on;
while BUDGET_MULTIPLIER*N > cocoProblemGetEvaluations(bbob_problem)
    
    feval_remain = BUDGET_MULTIPLIER*N-cocoProblemGetEvaluations(bbob_problem);
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

    %% APOP
    % [xmin,fmin,~] = apop_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, 1, SHOW_PLOW);

    %% LRA
    % [xmin,fmin,~] = lra_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, SHOW_PLOW);
    
    %% PSA
    % [xmin,fmin,~] = psa_cma_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, SHOW_PLOW);

    %% BIPOP
    % [xmin,fmin,~] = bipop_cma_es(fun, bbob_problem, mu, lambda, x0, lb, ub, sigma0, STOP, SHOW_PLOW);

    %% CMSA
    % [xmin,fmin,~] = cmsa_es(fun, bbob_problem, mu, lambda, x0, sigma0, STOP, SHOW_PLOW);

    %% META-ES
    lambda = 2^5;
    mu = lambda/2;
    [x,f,~] = meta_es(fun,bbob_problem, mu, lambda, x0,sigma0, STOP, SHOW_PLOW);

    if i<=NUM_PLOTS_SHOW
        fig = gcf;
        lineObjs = findobj(fig, 'type', 'line');
        if strcmpi(lineObjs(end).DisplayName, 'f')
            lineObjs(end).YData = lineObjs(end).YData-F_STOP;
        end
    end

    if cocoProblemFinalTargetHit(bbob_problem)
        break;
    end
    % plot(f_g-F_STOP, '-'); 
    % plot(sigD_g, 'k:'); 

    i = i+1;
end
%yscale('log');
%xlabel('g'); ylabel('f(g)-F');


if cocoProblemFinalTargetHit(bbob_problem)
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
