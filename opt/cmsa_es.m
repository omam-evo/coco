%% CMSA to call externally
% C = 1/sqrt(N_DIM);
% D = sqrt(N_DIM);
% E_chi = sqrt(N_DIM)*(1 - 1/(4*N_DIM) + 1/(21*N_DIM^2));

function [ymin, fmin, counteval, f_g, sig_g] = cmsa_es(fun,bbob,mu,lam,y,sigma,STOP,PRINT_SHOW)

    N = length(y);

    SET_WEIGHTED_REC = 1; 
    TAU = 1/sqrt(2*N); 
    TAU_C_MIN = 1;   % e.g. 4*N; needed for stability when MU/N is large

    ymin = y;
    fmin = fun(y);

    counteval = 1;  % '1' from fmin above
    g = 1;
    g_save = 1000;                 % dynamically append g_save values
    f_g = nan*zeros(g_save,1);     % dynamically append g_save values
    sig_g = nan*zeros(g_save,1);  % ...
    Dmax_g = nan*zeros(g_save,1);  % ...
    Dmin_g = nan*zeros(g_save,1);  % ...
    f_stag_g = nan*zeros(10,1);     % for stagnation check

    C = eye(N);
    sigmaTilde = nan*ones(1,lam);  dTilde = nan*ones(N,lam); 
    yTilde = nan*ones(N,lam); fTilde = nan*ones(1,lam);
    while counteval < STOP.FEVAL_MAX

        if SET_WEIGHTED_REC==1
            weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        else
            weights = 1/mu*ones(mu,1);
        end    
        weights = weights/sum(weights);         % normalize recombination weights array
        mueff = sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i
        TAU_C = max(1 + N*(N+1)/(2*mueff), TAU_C_MIN);
        y_old = y;

        %% EIG
        [V,D2] = eig(C);                             %eig()
        M = zeros(N,N);                             %eig()
        for i=1:N                                   %eig()
            M = M + sqrt(D2(i,i))*(V(:,i)*V(:,i)');  %eig()
        end                                         %eig()
        [min_D, max_D] = bounds(sqrt(diag(D2)));
        
        for l=1:lam
            sigmaTilde(l) = sigma * exp(TAU*randn());
            dTilde(:,l) = M*randn(N,1);  
            yTilde(:,l) = y + sigmaTilde(l)*dTilde(:,l);
            fTilde(l) = fun(yTilde(:,l));
        end
        counteval = counteval+lam;
        % -------------------------------------
        
        [fSorted, idx] =  sort(fTilde, 'ascend');
        y = yTilde(:, idx(1:mu)) * weights;
        sigma = sigmaTilde(idx(1:mu)) * weights;
        SumSS = zeros(N, N);
        for m=1:mu 
            SumSS = SumSS + dTilde(:, idx(m))*dTilde(:, idx(m))';
        end
        C = (1-1.0/TAU_C) * C + (1.0/TAU_C) * (1.0/mu) * SumSS;
        C = 0.5 * (C+C');
        fmin = fSorted(1);
        ymin = yTilde(:, idx(1));

        % save for stagnation check
        f_stag_g = stag_check(f_stag_g,fSorted,mu,g);

        %% Terminate
        if terminate_cma(bbob,STOP,fmin,g,sigma,max_D,min_D,y,y_old,f_stag_g)
            break;
        end

        %% Print results 
        if PRINT_SHOW==1
            disp(['g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fbest:' num2str(fSorted(1)), ' sig*D:' ... 
                  num2str(sigma*max_D,'%.4e'), ' cond-sqrt:' ...
                  num2str(max_D/min_D,'%.2e'), ' mu:', num2str(mu,'%i')]);
            f_g(g) = fSorted(1);
            sig_g(g) = sigma*max_D;
            Dmax_g(g) = max_D;
            Dmin_g(g) = min_D;
            if mod(g, g_save)==0 % preallocate g_save slots
                f_g = [f_g; nan*zeros(g_save,1)];
                sig_g = [sig_g; nan*zeros(g_save,1)];
                Dmax_g = [Dmax_g; nan*zeros(g_save,1)];
                Dmin_g = [Dmin_g; nan*zeros(g_save,1)];
            end
        end
        
        g = g+1;

    end % counteval
    if counteval >= STOP.FEVAL_MAX
        fprintf('\t FEVAL STOP STOP.\n');
    end

    if PRINT_SHOW==1
        figure; hold on; title('cmsa-es');
            plot(f_g, 'DisplayName','f');
            plot(sig_g, 'DisplayName','sigma');       
            plot(Dmax_g, 'DisplayName','sqrt D max');  
            plot(Dmin_g, 'DisplayName','sqrt D min');  
            yscale('log'); legend;
    end

end
