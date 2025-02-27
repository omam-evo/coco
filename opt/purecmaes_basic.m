function [xmin,fmin,counteval,xmean,sigma,C,max_D,min_D] = purecmaes_basic(fun, bbob, mu, lambda, xmean, sigma, STOP, PRINT_SHOW)

% (mu/mu_w, lambda)-CMA-ES 
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation
% for nonlinear function minimization. 
%
% This code is "an excerpt" from cmaes.m and implements the key 
% parts of the algorithm. It is intendend to be used for READING
% and UNDERSTANDING the basic flow and all details of the CMA-ES
% *algorithm*. To run "serious" simulations better use the cmaes.m 
% code: it is longer, but offers restarts, far better termination 
% options, and, in particular, supposedly quite useful output.
%
% Author: Nikolaus Hansen, 2003-09. 
% e-mail: hansen[at]lri.fr
%
% License: This code is released into the public domain (that is, 
%   you may use and modify it however you like). 
%
% URL: http://www.lri.fr/~hansen/purecmaes.m
% References: See end of file. Last change: April, 29, 2014

% --------------------  Initialization --------------------------------  
% User defined input parameters (need to be edited)

    CMA_ON = 1;           % default: 1  (only CSA: 0) 
    CSA_MODE = 1;         % default: 1  (alternative: 2)
    SET_WEIGHTED_REC = 1; % default: 1  (intermediate: 0) 
    
    N = length(xmean);

    xmin = xmean;
    fmin = fun(xmin);
    
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
    sigma_D_g = nan*zeros(g_save,1);
    Dmax_g = nan*zeros(g_save,1);
    Dmin_g = nan*zeros(g_save,1);
    f_stag_g = nan*zeros(10,1);

    while counteval < STOP.FEVAL_MAX
        
        % Generate and evaluate lambda offspring
        for k=1:lambda
          arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
          arfitness(k) = fun(arx(:,k)); % objective function call
        end
        counteval = counteval+lambda;
        
        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness);  % minimization
        xold = xmean;
        xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
        xmin = arx(:, arindex(1));
        fmin = arfitness(1);
        
        % save for stagnation check
        f_stag_g = stag_check(f_stag_g,arfitness,mu,g);
        
        % Cumulation: Update evolution paths
        ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
 
        % Adapt covariance matrix C
        if CMA_ON==1 % && lambda > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
            %hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
            hsig = 1;
            pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma; 
            artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
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
        if terminate_cma(bbob,STOP,fmin,g,sigma,max_D,min_D,xmean,xold,f_stag_g)
            break;
        end
    
        %% Output 
        if PRINT_SHOW==1
            disp([num2str(counteval,'%i') ': ' num2str(fmin) ' ' ... 
                  num2str(sigma*max_D,'%.4e') ' ' ...
                  num2str(max_D/min_D,'%.2e')]);
            f_g(g) = fmin;
            sigma_D_g(g) = sigma*max_D;
            Dmax_g(g) = max_D;
            Dmin_g(g) = min_D;
            if mod(g, g_save)==0 % preallocate g_save slots
                f_g = [f_g; nan*zeros(g_save,1)];
                Dmax_g = [Dmax_g; nan*zeros(g_save,1)];
                Dmin_g = [Dmin_g; nan*zeros(g_save,1)];
                sigma_D_g = [sigma_D_g; nan*zeros(g_save,1)];
            end
        end
        g = g+1;
    
        % with long runs, the next line becomes time consuming
        % out.dat = [out.dat; arfitness(1) sigma 1e5*D' ]; 
        % out.datx = [out.datx; xmean'];

    end % while, end generation loop


    if PRINT_SHOW==1
        figure; hold on; title('purecmaes');
            plot(f_g, 'DisplayName','f');
            plot(sigma_D_g, 'DisplayName','sigma');   
            plot(Dmax_g, '--', 'DisplayName','Dmax');   
            plot(Dmin_g, '--', 'DisplayName','Dmin');   
            set(gca,'Yscale','log'); legend;
    end

end % function
