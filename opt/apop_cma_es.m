function [xmin,fmin,f_g,sigD_g,lam_g] = apop_cma_es(fun, bbob, mu, lambda, x_mean, sigma, STOP, APOP_ON, PRINT_SHOW)

    CMA_ON = 1;           % default: 1  (only CSA: 0) 
    CSA_MODE = 1;         % default: 1  (alternative: 2)
    SET_WEIGHTED_REC = 1; % default: 1  (intermediate: 0) 
    % Learning Rates by Nomura & Akimoto: eta_m,eta_s            

    N = length(x_mean);
    xmin = x_mean;
    fmin = fun(xmin);
    
    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
    B = eye(N,N);                       % B defines the coordinate system
    D = ones(N,1);                      % diagonal D defines the scaling
    C = B * diag(D.^2) * B';            % covariance matrix C
    invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
    chiN = N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of ||N(0,I)|| == norm(randn(N,1))
    % eigeneval = 0;                      % track update of B and D
    % out.dat = []; out.datx = [];  % for plotting output
    
    % -------------------- Generation Loop --------------------------------
    counteval = 1;  % '1' from fmin above; the next 40 lines contain the 20 lines of interesting code 
    g = 1;
    g_save = 1000;
    f_g = nan*zeros(g_save,1);
    sigma_g = nan*zeros(g_save,1);
    sigD_g = nan*zeros(g_save,1);
    lam_g = nan*zeros(g_save,1);

    %% APOP
    THETA = 1/2;
    lam_def = 4+floor(3*log(N));  
    MU_MIN = 2*lam_def*THETA;   % Line 30 [NH17]
    MU_MAX = 400*lam_def*THETA; % Line 25 [NH17]
    L = 6;                      % L=6 yields L-1=5=S f-differences, Lines 7-13, S=5 below Algorithm
    P_Thr = 0.2;                % 1/5, Line 30 [NH17]
    wait = L;  
    alpha_mu = 2;
    f_med_g = nan*ones(1,L);
    f_stag_g = nan*zeros(10,1);

    while counteval < STOP.FEVAL_MAX
        
        % Strategy parameter setting: Selection  
        if SET_WEIGHTED_REC==1
            weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        else
            weights = 1/mu*ones(mu,1);
        end    
        weights = weights/sum(weights);         % normalize recombination weights array
        mueff = sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i
        
        % Strategy parameter setting: Adaptation
        if CMA_ON==1
            cc = (4 + mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
            c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
            cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
        else
            c1 = 0; cc = 0; cmu = 0;
        end
        if CSA_MODE == 1
            cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
            damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs;        % damping for sigma 
        elseif CSA_MODE == 2
            cs = 1/sqrt(N);
            damps = 1/cs;
        end

        % Generate and evaluate lambda offspring
        for k=1:lambda
          arx(:,k) = x_mean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
          arfitness(k) = fun(arx(:,k)); % objective function call
        end
        counteval = counteval+lambda;
        
        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness);  % minimization
        x_old = x_mean;
        xmin = arx(:, arindex(1));
        fmin = arfitness(1);

        % save for stagnation check
        f_stag_g = stag_check(f_stag_g,arfitness,mu,g);
    
        % Update mean vector
	    x_mean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
    
        % Cumulation: Update evolution paths
        ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (x_mean-x_old) / sigma; 
 
        % Adapt covariance matrix C
        if CMA_ON==1 % && lambda > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
            %hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
            hsig = 1;
            pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (x_mean-x_old) / sigma; 
            artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(x_old,1,mu));  % mu difference vectors
            C = (1-c1-cmu) * C ...                   % regard old matrix  
                 + c1 * (pc * pc' ...                % plus rank one update
                         + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
                 + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 
        end
    
        % Adapt step size sigma
        if CSA_MODE == 1
            sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
        elseif CSA_MODE == 2
            sigma = sigma * exp((norm(ps)/chiN - 1)/damps); 
        end
        % Update B and D from C
        if CMA_ON==1 %&& lambda > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
          C = triu(C) + triu(C,1)'; % enforce symmetry
          [B,D2] = eig(C);           % eigen decomposition, B==normalized eigenvectors
	      D = sqrt(diag(D2));        % D contains standard deviations now
          invsqrtC = B * diag(D.^-1) * B';
        end
        [min_D,max_D] = bounds(D);

        %% Terminate
        if terminate_cma(bbob,STOP,fmin,g,sigma,max_D,min_D,x_mean,x_old,f_stag_g)
            break;
        end
    
        %% Output before population is changed
        if PRINT_SHOW==1
            disp(['g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fmin:' num2str(fmin), ' sig*D:' ... 
                  num2str(sigma*max_D,'%.4e'), ' cond:' ...
                  num2str(max_D/min_D,'%.2e'), ' lambda:', num2str(lambda,'%i')]);
            f_g(g) = arfitness(1);
            sigD_g(g) = sigma*max_D;
            sigma_g(g) = sigma;
            lam_g(g) = lambda;
            if mod(g, g_save)==0 % preallocate g_save slots
                f_g = [f_g; nan*zeros(g_save,1)];
                sigD_g = [sigD_g; nan*zeros(g_save,1)];
                sigma_g = [sigma_g; nan*zeros(g_save,1)];
            end
        end

        %% Population Size Control
        if APOP_ON == 1 
            if g>=L   % g>=L: APOP evaluates fitness differences
                df = diff(f_med_g);
                P = sum(df>0)/(L-1);
            end
            mu_old = mu;
            if wait <= 0
                % if P==P_Thr     % neutral performance
                if P<P_Thr  % positive performance
                    mu = mu/alpha_mu;
                    mu = floor(mu);
                elseif P>P_Thr  % negative performance
                    mu = mu*alpha_mu;
                    mu = ceil(mu);
                end
                mu = min(max(mu, MU_MIN), MU_MAX); % keep in bounds
            end 
           
            % new mu, adjust params
            if mu~=mu_old 
                wait = L; % wait only if pop. change happened
                sigma = sigma*(mu/mu_old)^0.5;  
                lambda = floor(mu*2);
            else % no change, decrease wait
                wait = wait-1;
            end
           
            if g<=L
                f_med_g(g) = median(arfitness(1:mu_old));
            else
                f_med_g = [f_med_g(2:end), median(arfitness(1:mu_old))];
            end
        end % APOP_ON

        g = g+1;
    
        % with long runs, the next line becomes time consuming
        % out.dat = [out.dat; arfitness(1) sigma 1e5*D' ]; 
        % out.datx = [out.datx; xmean'];

    end % while, end generation loop

    if PRINT_SHOW==1
        figure; hold on;
            plot(f_g, 'DisplayName','f');
            plot(lam_g, 'DisplayName','lambda');
            %plot(sigD_g, 'DisplayName','sigD');  
            plot(sigma_g, 'DisplayName','sigma'); 
            yscale('log'); legend;
    end

end % function
