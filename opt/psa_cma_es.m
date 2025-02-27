%% Implementation of PSA-CMA-ES from Nishida, Akimoto 2018 (CMA from Hansen purecmaes.m)
% the following simplifications were implemented since the original algorithm is involved
% > "hsig=1" was set (in default CMA-ES hsig=0 if length of cumulation path becomes very large, pausing update of pc)
% > SET_gamma_time_dep=0 simplifies time-dependent gamma (normalization) quantities and E_I
% (reason: original PSA includes time-dependent normalization factors gamma_sigma, gamma_c, gamma_theta, of which the effects should be "barely recognizable in practice", page 3, top)
% > simplified sigma-rescaling using ratio of mu_w (see (18) in PSA paper assuming mu_w*n>>mu_w); hence, coeff. c not needed
% > E_I from (13) was also simplified including only the first term n/mu_w
% > default CMA-ES uses many more termination criteria, which were not implemented here

function [xmin,fmin,counteval,xmean,sigma,C,max_D,min_D] = psa_cma_es(fun, bbob, mu, lambda, xmean, sigma, STOP, PRINT_SHOW)

    CMA_ON = 1;           % default: 1  (only CSA: 0) 
    CSA_MODE = 1;         % default: 1  (alternative: 2)
    SET_WEIGHTED_REC = 1; % default: 1  (intermediate: 0) 

    N = length(xmean);

    xmin = xmean;
    fmin = fun(xmin);
    
    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
    B = eye(N,N);                       % B defines the coordinate system
    D = ones(N,1);                      % diagonal D defines the scaling
    C = B * diag(D.^2) * B';            % covariance matrix C
    invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
    chiN = N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of ||N(0,I)|| == norm(randn(N,1))

    % -------------------- Generation Loop --------------------------------
    counteval = 1;  % '1' from fmin above; the next 40 lines contain the 20 lines of interesting code 
    g = 1;
    g_save = 1000;
    f_g = nan*zeros(g_save,1);
    sigD_g = nan*zeros(g_save,1);
    sig_g = nan*zeros(g_save,1);
    Dmax_g = nan*zeros(g_save,1);
    Dmin_g = nan*zeros(g_save,1);
    pm_g = nan*zeros(g_save,1);
    ps_g = nan*zeros(g_save,1);
    lam_g = nan*zeros(g_save,1);
    f_stag_g = nan*zeros(10,1); 

    alpha = 1.4;
    beta = 0.4;
    SET_gamma_time_dep = 0;
    if SET_gamma_time_dep == 1
        gamma_c = 0; 
        gamma_s = 0;
        gamma_t = 0; 
    else
        gamma_c = 1; 
        gamma_s = 1;
        gamma_t = 1;    
    end
    psa_m = zeros(N,1);
    psa_s = zeros(N+N*(N-1)/2,1);
    lambda_max = 512*(4+floor(3*log(N)));
    lambda_min = lambda;

    while counteval < STOP.FEVAL_MAX

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
            ds = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs;        % damping for sigma 
        elseif CSA_MODE == 2
            cs = 1/sqrt(N);
            ds = 1/cs;
        end

        xmean_old = xmean;
        sigma_old = sigma;
        C_old = C;
        lambda_old = lambda;

        % Generate and evaluate lambda offspring
        for k=1:lambda
          arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
          arfitness(k) = fun(arx(:,k)); % objective function call
        end
        counteval = counteval+lambda;
        
        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness);  % minimization
        xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
        xmin = arx(:, arindex(1));
        fmin = arfitness(1);

        % save for stagnation check
        f_stag_g = stag_check(f_stag_g,arfitness,mu,g);
    
        % Cumulation: Update evolution paths
        ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xmean_old) / sigma; 
 
        % Adapt covariance matrix C

        if CMA_ON==1
            % hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
            hsig = 1;
            pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xmean_old) / sigma; 
            artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xmean_old,1,mu));  % mu difference vectors
            
            %PSA
            if SET_gamma_time_dep == 1
                gamma_s = (1-cs)^2*gamma_s + cs*(2-cs);
                gamma_c = (1-cc)^2*gamma_c + hsig*cc*(2-cc);
            end
            
            C = (1-c1*gamma_c-cmu) * C ...                   % regard old matrix  
                 + c1 * (pc * pc' ...                % plus rank one update
                         + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
                 + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 
        end

        if SET_gamma_time_dep == 1
            E_I = N/mueff + 2*N*(N-chiN^2)/chiN^2*gamma_s*(cs/ds)^2;
            E_I = E_I + 0.5*(1+8*gamma_s*(N-chiN^2)/chiN^2*(cs/ds)^2)* ...
                ((N^2+N)*cmu^2/mueff + (N^2+N)*cc*(2-cc)*c1*cmu*mueff*sum(weights.^3)+ ...
                c1^2* (gamma_c^2*N^2 + N*(1-2*gamma_c+2*gamma_c^2) ) );
        else
            E_I = N/mueff;
        end
    
        % Adapt step size sigma
        if CSA_MODE == 1
            sigma = sigma * exp((cs/ds)*(norm(ps)/chiN - gamma_s)); 
        elseif CSA_MODE == 2
            sigma = sigma * exp((norm(ps)/chiN - 1)/ds); 
        end
        % Update B and D from C
        if CMA_ON==1 %&& lambda > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
          C = triu(C) + triu(C,1)'; % enforce symmetry
          [B,D2] = eig(C);           % eigen decomposition, B==normalized eigenvectors
	      D = sqrt(diag(D2));        % D contains standard deviations now
          invsqrtC = B * diag(D.^-1) * B';
        end
        [min_D,max_D] = bounds(D);

        %% PSA
        dm = xmean - xmean_old;
        S_old = sigma_old^2*C_old;
        S = sigma^2*C;
        dS = S - S_old;

        % normalize
        invsqrtS = (1/sigma)*invsqrtC;
        dm_tilde = invsqrtS * dm;  
        dS_mat = 2^(-1/2)*(invsqrtS*dS*invsqrtS);
        dS_tilde = dS_mat(triu(dS_mat)~=0);

        psa_m = (1-beta)*psa_m + sqrt(beta*(2-beta))*dm_tilde/sqrt(E_I);
        psa_s = (1-beta)*psa_s + sqrt(beta*(2-beta))*dS_tilde/sqrt(E_I);
        psa_m2 = psa_m'*psa_m;
        psa_s2 = psa_s'*psa_s;
        psa_t2 = psa_m2 + psa_s2;
        if SET_gamma_time_dep == 1
            gamma_t = (1-beta)^2*gamma_t + beta*(2-beta);
        end

        lambda = lambda*exp(beta*(gamma_t-psa_t2/alpha));
        lambda = round(min(max(lambda,lambda_min),lambda_max));
        mu = floor(lambda/2);

        %% Linear sigma-rescaling
        if SET_WEIGHTED_REC==1
            weights_new = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        else
            weights_new = 1/mu*ones(mu,1);
        end    
        weights_new = weights_new/sum(weights_new);         % normalize recombination weights array
        mueff_new = sum(weights_new)^2/sum(weights_new.^2);
        sigma = sigma*(mueff_new/mueff);
        
        %% Terminate
        if terminate_cma(bbob,STOP,fmin,g,sigma,max_D,min_D,xmean,xmean_old,f_stag_g)
            break;
        end
    
        %% Output 
        if PRINT_SHOW==1
            disp(['g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fbest:' num2str(fmin), ' sig*D:' ... 
                  num2str(sigma*max_D,'%.4e'), ' cond-sqrt:' ...
                  num2str(max_D/min_D,'%.2e'), ' pm:', num2str(psa_m2,'%.4e'), ' ps:', num2str(psa_s2,'%.4e')]);
            f_g(g) = arfitness(1);
            sig_g(g) = sigma;
            Dmax_g(g) = max_D;
            Dmin_g(g) = min_D;
            pm_g(g) = psa_m2;
            ps_g(g) = psa_s2;
            lam_g(g) = lambda;
            if mod(g, g_save)==0 % preallocate g_save slots
                f_g = [f_g; nan*zeros(g_save,1)];
                Dmax_g = [Dmax_g; nan*zeros(g_save,1)];
                Dmin_g = [Dmin_g; nan*zeros(g_save,1)];
                sig_g = [sig_g; nan*zeros(g_save,1)];
                pm_g = [pm_g; nan*zeros(g_save,1)];
                ps_g = [ps_g; nan*zeros(g_save,1)];
                lam_g = [lam_g; nan*zeros(g_save,1)];
            end
        end
        g = g+1;
    
        % with long runs, the next line becomes time consuming
        % out.dat = [out.dat; arfitness(1) sigma 1e5*D' ]; 
        % out.datx = [out.datx; xmean'];

    end % while, end generation loop


    if PRINT_SHOW==1
        figure; hold on; title('psa-cma-es');
            plot(f_g, 'DisplayName','f');
            plot(sig_g, 'DisplayName','sigma');       
            plot(pm_g+ps_g, 'DisplayName','pt');
            %plot(ps_g, 'DisplayName','ps2');
            plot(lam_g, 'k-', 'DisplayName','lambda');
            yscale('log'); legend;
    end

end % function
