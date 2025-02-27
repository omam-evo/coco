function [x,f,counteval,f_g,lambda_g] = meta_es(fun,bbob,mu,lambda,ISO_FACTOR,x,sigma,STOP,PRINT_SHOW)
    
    PRINT_SHOW_INNER = 0;
    THETA = mu/lambda;
    N = length(x);

    %% Population
    POP_FACTOR = 2;
    LAMBDA_MIN = 8; %max(8,2^floor(log2(N)));
    LAMBDA_MAX = 4096;
    MODE_ISOLATION = 1; %{1,2,3}
    
    %% Init
    C = eye(N);         % CMA-ES would require: ps = zeros(N,1); pc = zeros(N,1); 
    STOP_INNER = STOP;  % inner stopping critera are modified below
    counteval = 0;
    g = 1;
    f_g = [];       % dynamically append f of sub-results
    lambda_g = [];  % dynamically append pop. size

    %ISO_FACTOR = 10;   % isolation factor (lambda2)
    f_stag_g = nan*zeros(10,1); 

    while counteval < STOP.FEVAL_MAX

        x_old = x;
        lambda1 = max(ceil(lambda/POP_FACTOR), LAMBDA_MIN);
        lambda2 = min(ceil(lambda*POP_FACTOR), LAMBDA_MAX);
        r_sigma = sqrt(POP_FACTOR);

        if MODE_ISOLATION==1
            STOP_INNER.FEVAL_MAX = ISO_FACTOR*lambda2; % dynamic func.eval: ensure larger pop. runs for multiple generations
        elseif MODE_ISOLATION==2
            STOP_INNER.FEVAL_MAX = LAMBDA_MAX; % fixed func.eval: max. pop. will run only for 1 generation
        elseif MODE_ISOLATION==3
            STOP_INNER.G_MAX = ISO_FACTOR; % fix number gen. for isolation => pop. size goes up on unimodal func.
        end
        %% Remarks:
        % inner ES have only function evaluations as termination condition
        % inner ES uses lambda+1 function evaluations per generation (+1 for recombinant fitness evaluation: improve global search)
        % ISO_FACTOR*lambda2 ensures that larger population has certain number of gen., while smaller has roughly 4 times the budget
        [x1, f1, counteval1, sigma1, C1, max_D1, min_D1, f1_g] = ...
            cmsa_es_meta(fun,bbob,ceil(lambda1*THETA),lambda1,x,sigma/r_sigma,C,STOP_INNER,PRINT_SHOW_INNER);
        [x2, f2, counteval2, sigma2, C2, max_D2, min_D2, f2_g] = ...
            cmsa_es_meta(fun,bbob,ceil(lambda2*THETA),lambda2,x,sigma*r_sigma,C,STOP_INNER,PRINT_SHOW_INNER);

        counteval = counteval + counteval1 + counteval2;

        %% Display fitness dynamics as function of function evaluations
        % gen1 = sum(~isnan(f1_g));
        % gen2 = sum(~isnan(f2_g));h
        %     plot((g-1)*isolate_feval+mu1/THETA*(0:gen1-1),f1_g(1:gen1), 'b-', 'DisplayName',['mu1=',num2str(mu1)]);
        %     plot((g-1)*isolate_feval+mu2/THETA*(0:gen2-1),f2_g(1:gen2), 'r-.', 'DisplayName',['mu2=',num2str(mu2)]);
        %     xline(g*isolate_feval, ':')
        %     yscale('log'); %legend('Location','northoutside');

        if f1<=f2
            x = x1; f=f1; sigma = sigma1; lambda=lambda1; C = C1; max_D = max_D1; min_D = min_D1;
            f_g = [f_g;f1_g];
            lambda_g = [lambda_g;lambda*ones(length(f1_g),1)];
        else
            x = x2; f=f2; sigma = sigma2; lambda=lambda2; C = C2; max_D = max_D2; min_D = min_D2;
            f_g = [f_g;f2_g];
            lambda_g = [lambda_g;lambda*ones(length(f2_g),1)];
        end

        % save for stagnation check
        f_stag_g = stag_check(f_stag_g,f,1,g);

        if PRINT_SHOW==1
            disp(['meta-es: g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fbest:' num2str(f), ' sig*D:' ... 
                  num2str(sigma*max_D,'%.4e'), ' cond-sqrt:' ...
                  num2str(max_D/min_D,'%.2e'), ' lambda:', num2str(lambda,'%i')]);
        end

        %% Terminate
        if terminate_cma(bbob,STOP,f,g,sigma,max_D,min_D,x,x_old,f_stag_g)
            break;
        end

        g = g+1;
    end
    if PRINT_SHOW==1
        figure; hold on; title('meta-es');
            plot(f_g, 'DisplayName','f');
            plot(lambda_g, 'DisplayName','lambda');  
            yscale('log'); legend;
    end

end