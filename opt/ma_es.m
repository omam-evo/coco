%% CMSA to call externally
% C = 1/sqrt(N_DIM);
% D = sqrt(N_DIM);
% E_chi = sqrt(N_DIM)*(1 - 1/(4*N_DIM) + 1/(21*N_DIM^2));

function [ymin, fmin, counteval, f_g, sig_g] = ma_es(fun,bbob,mu,lam,y,sigma,STOP,PRINT_SHOW)

    N = length(y);

    SET_WEIGHTED_REC = 1; 

    ymin = y;
    fmin = fun(y);

    counteval = 1;  % '1' from fmin above
    g = 1;
    g_save = 100;                 % dynamically append g_save values
    f_g = nan*zeros(g_save,1);     % dynamically append g_save values
    sig_g = nan*zeros(g_save,1);  % ...

    M = eye(N);
    I = eye(N);
    chiN = N^0.5*(1-1/(4*N)+1/(21*N^2)); 
    s = zeros(N,1);
    zTilde = nan*ones(N,lam);  dTilde = nan*ones(N,lam); 
    yTilde = nan*ones(N,lam); fTilde = nan*ones(1,lam);

    if SET_WEIGHTED_REC==1
        weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
    else
        weights = 1/mu*ones(mu,1);
    end    
    weights = weights/sum(weights);         % normalize recombination weights array
    mueff = sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

    c1 = 2 / ((N+1.3)^2+mueff); 
    cc = (4 + mueff/N) / (N+4 + 2*mueff/N);
    cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));
    cs = (mueff+2) / (N+mueff+5); 
    ds = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; 

    while counteval < STOP.FEVAL_MAX

        %% EIG
        %[min_D, max_D] = bounds(sqrt(diag(D2)));
        
        for l=1:lam
            zTilde(:,l) = randn(N,1);  
            dTilde(:,l) = M*zTilde(:,l);
            yTilde(:,l) = y + sigma*dTilde(:,l);
            fTilde(l) = fun(yTilde(:,l));
        end
        counteval = counteval+lam;
        
        [fSorted, idx] =  sort(fTilde, 'ascend');
        fmin = fSorted(1);
        ymin = yTilde(:, idx(1));

        y_old = y;
        y = y_old + sigma*(dTilde(:, idx(1:mu)) * weights);
        s = (1-cs) * s + sqrt(cs*(2-cs)*mueff) * (zTilde(:, idx(1:mu)) * weights); 
        
        zz = zeros(N,N);
        for m=1:mu
            zz = zz + weights(m)*zTilde(:, idx(m))*zTilde(:, idx(m))';
        end

        M = M*(I + 0.5*c1*(s*s'-I) + 0.5*cmu*(zz-I));

        sigma = sigma * exp((cs/ds)*(norm(s)/chiN - 1)); 

       %% Terminate
        if ~isempty(bbob)
            if cocoProblemFinalTargetHit(bbob) %cocoProblemGetBestValue(problem)
                fprintf('SUCCESS cocoProblemFinalTargetHit.\n');
                break
            end
        elseif fSorted(1) <= STOP.F_STOP % not needed for coco
            fprintf('\t F STOP.\n');
            break;
        end
        if norm(y-y_old) < STOP.X_TOL % sigma and maximum eigenvalue stop
            fprintf('\t SIGMA*D STOP.\n');
            break;
        % if sigma*max_D < SIGMA_D_STOP % sigma and maximum eigenvalue stop
        %     fprintf('\t SIGMA*D STOP purecmaes.\n');
        %     break;
        % elseif max_D/min_D > STOP.COND_SQRT %sqrt of condition 1e14 => 1e7
        %     fprintf('\t COND STOP.\n');
        %     break;
        end
        if PRINT_SHOW==1
            disp(['g:',num2str(g,'%i'), ' feval:', num2str(counteval,'%i'), ' fbest:' num2str(fSorted(1)), ' mu:', num2str(mu,'%i')]);
            f_g(g) = fSorted(1);
            sig_g(g) = sigma;
            %Dmax_g(g) = max_D;
            %Dmin_g(g) = min_D;
            if mod(g, g_save)==0 % preallocate g_save slots
                f_g = [f_g; nan*zeros(g_save,1)];
                sig_g = [sig_g; nan*zeros(g_save,1)];

                [D2] = eig(M*M');
                [min_D,max_D] = bounds(sqrt(D2));
                Dmax_g(g-g_save+1:g,1) = max_D;
                Dmin_g(g-g_save+1:g,1) = min_D;
            end
        end
        
        g = g+1;

    end % counteval
    if counteval >= STOP.FEVAL_MAX
        fprintf('\t FEVAL STOP STOP.\n');
    end

    if PRINT_SHOW==1
        figure; hold on; title('ma-es');
            plot(f_g, 'DisplayName','f');
            plot(sig_g, 'DisplayName','sigma'); 
            plot(Dmax_g, '--', 'DisplayName','sqrt D max');  
            plot(Dmin_g, '--', 'DisplayName','sqrt D min');  
            yscale('log'); legend;
    end

end
