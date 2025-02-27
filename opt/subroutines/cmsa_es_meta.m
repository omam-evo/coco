%% CMSA to call externally
% C = 1/sqrt(N_DIM);
% D = sqrt(N_DIM);
% E_chi = sqrt(N_DIM)*(1 - 1/(4*N_DIM) + 1/(21*N_DIM^2));

function [y, f, counteval, sigma, C, max_D, min_D, f_g] = cmsa_es_meta(fun,bbob,mu,lam,y,sigma,C,STOP,PRINT_SHOW)

    N = length(y);

    SET_WEIGHTED_REC = 1; 
    TAU = 1/sqrt(2*N);   %TAU = 1/sqrt(2*N); 
    TAU_C_MIN = 1; % e.g. 4*N; may be needed for stability when MU/N is large

    f = nan;

    counteval = 0;  % '1' from fmin above
    g = 1;
    g_save = 1000;                 % dynamically append g_save values
    f_g = nan*zeros(g_save,1);     % dynamically append g_save values
    sigma_g = nan*zeros(g_save,1);  % ...
    Dmax_g = nan*zeros(g_save,1);  % ...
    Dmin_g = nan*zeros(g_save,1);  % ...
    % fmed_g = nan*zeros(10,1);   % stagnation check

    %C = eye(N); % input
    sigmaTilde = nan*ones(1,lam);  dTilde = nan*ones(N,lam); 
    yTilde = nan*ones(N,lam); fTilde = nan*ones(1,lam);
    while counteval+lam <= STOP.FEVAL_MAX && g <= STOP.G_MAX

        y_old = y;

        if SET_WEIGHTED_REC==1
            weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        else
            weights = 1/mu*ones(mu,1);
        end    
        weights = weights/sum(weights);         % normalize recombination weights array
        mueff = sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i
        TAU_C = max(1 + N*(N+1)/(2*mueff), TAU_C_MIN);

        %% EIG
        [V,D2] = eig(C);                            %eig()
        M = zeros(N,N);                             %eig()
        for i=1:N                                   %eig()
            M = M + sqrt(D2(i,i))*(V(:,i)*V(:,i)'); %eig()
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

        %% Choose best or recombinant:
        % f = fSorted(1);
        % y = yTilde(:, idx(1));
        f = fun(y); counteval = counteval+1;

        % Check stagnation
        % if g>10
        %     fmed_g(1:end-1) = fmed_g(2:end);
        %     fmed_g(end) = median(fSorted);
        % else
        %     fmed_g(g) = median(fSorted);
        % end


       %% Terminate
        % if ~isempty(bbob)
        %     if cocoProblemFinalTargetHit(bbob) %cocoProblemGetBestValue(problem)
        %         fprintf('SUCCESS cocoProblemFinalTargetHit.\n');
        %         break
        %     end
        % elseif fSorted(1) <= STOP.F_STOP % not needed for coco
        %     fprintf('   F STOP.\n');
        %     break;
        % end
        % if sigma*max_D < STOP.SIGMA_D_STOP % sigma and maximum eigenvalue stop
        %     fprintf('   SIGMA*D STOP cmsa.\n');
        %     break;
        % elseif max_D/min_D > STOP.COND_SQRT %sqrt of condition 1e14 => 1e7
        %     fprintf('   COND STOP cmsa.\n');
        %     break;
        % elseif norm(y-y_old) < STOP.X_TOL
        %     fprintf('   XTOL STOP cmsa.\n');
        %     break;   
        % elseif g>10 && max(fmed_g)-min(fmed_g) < STOP.F_STAG_TOL
        %     fprintf('   F_STAG_TOL STOP cmsa.\n');
        %     break;  
        % end

        if mod(g, g_save)==0 % preallocate g_save slots
            f_g = [f_g; nan*zeros(g_save,1)];
            sigma_g = [sigma_g; nan*zeros(g_save,1)];
            Dmax_g = [Dmax_g; nan*zeros(g_save,1)];
            Dmin_g = [Dmin_g; nan*zeros(g_save,1)];
        end
        f_g(g) = f;
        sigma_g(g) = sigma*max_D;
        Dmax_g(g) = max_D;
        Dmin_g(g) = min_D;

        if PRINT_SHOW==1
            disp(['   g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fbest:' num2str(fSorted(1)), ' sig*D:' ... 
                  num2str(sigma*max_D,'%.4e'), ' cond-sqrt:' ...
                  num2str(max_D/min_D,'%.2e'), ' lambda:', num2str(lam,'%i')]); 
        end
        
        g = g+1;

    end % counteval
    f_g = f_g(~isnan(f_g));

    if PRINT_SHOW==1 && counteval+lam > STOP.FEVAL_MAX
            disp(['   Isolation:', ' feval:', num2str(counteval,'%i'), ' f:' num2str(f), ' lambda:', num2str(lam,'%i')]);
    end

    if PRINT_SHOW==1
        figure; hold on; title('cmsa-es');
            plot(f_g, 'DisplayName','f');
            plot(sigma_g, 'DisplayName','sigma');       
            plot(Dmax_g, 'DisplayName','sqrt D max');  
            plot(Dmin_g, 'DisplayName','sqrt D min');  
            yscale('log'); legend;
    end

end
