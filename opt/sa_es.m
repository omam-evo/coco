function [y, f, sigma, counteval, f_g, sigma_g] = ...
    sa_es(fun, bbob, mu, lam, y, sigma, TAU, STOP, PRINT_SHOW)

    N = length(y);

    % Values over generations
    g = 1;
    g_save = 1000;
    f_g = nan*zeros(g_save,1); 
    sigma_g = nan*ones(g_save,1); 
    r_g = nan*ones(g_save,1); 
    y_tilde = nan*ones(N, lam); 
    s_tilde = nan*ones(1, lam); 
    f_tilde = nan*ones(1, lam);
    counteval = 0;
    
    while counteval < STOP.FEVAL_MAX 

        %% Sequential

        for L = 1:lam
            s_tilde(1,L) = sigma * exp(TAU * randn(1));    % LOGNORMAL
            % s_tilde(1,L) = sigma * (1 + TAU * randn(1)); % NORMAL
            y_tilde(:, L) = y + s_tilde(1,L) * randn(N,1);
            f_tilde(1,L) = fun(y_tilde(:, L));
        end

        %%   Vectorized
        % s_tilde = repmat(sigma * (1 + TAU * randn(1, LAM)), N, 1);    % NORMAL
        % % s_tilde = repmat(sigma * exp(TAU * randn(1, LAM)), N, 1); % LOGNORMAL
        % y_tilde = repmat(y, 1, LAM) + s_tilde.*randn(N, LAM);
        % [f_tilde] = fun(y_tilde); 
        % 
        counteval = counteval + lam;

        %% Selection
        [f_sorted, idx] = sort(f_tilde, 'ascend');
        y_n_mu = y_tilde(:, idx(1:mu));
        s_mu = s_tilde(1, idx(1:mu));

        %% Intermediate Recombination
        y = mean(y_n_mu, 2);
        sigma = mean(s_mu, 2);
        f = fun(y);
        fmin = f_sorted(1);
        counteval = counteval + 1;

        %% Save data
        if mod(g, g_save)==0 % pPRINT_GENreallocate g_save slots
            f_g = [f_g; nan*zeros(g_save,1)];
            r_g = [r_g; nan*zeros(g_save,1)];
            sigma_g = [sigma_g; nan*zeros(g_save,1)];
        end
        f_g(g) = f;
        sigma_g(g) = sigma;
        r_g(g) = norm(y);
        
        %% Print
        if PRINT_SHOW == 1
            fprintf('g: %i, feval: %i, f_rec: %.6e, sigma: %.6e,  lam: %i \n', g, counteval, f, sigma, lam);
        end
        %% Terminate
        if ~isempty(bbob)         % if script is called with bbob-problem
            if cocoProblemFinalTargetHit(bbob) %cocoProblemGetBestValue(problem)
                fprintf('SUCCESS.\n');
                break;
            end
        elseif f <= STOP.F_STOP % if script is called with F_STOP (no bbob)
            fprintf('F STOP.\n');
            break;
        end
        if sigma < STOP.SIGMA_STOP
            if PRINT_SHOW==1
                fprintf('SIGMA STOP.\n');
            end
            break
        elseif g == STOP.G_MAX
            fprintf('GEN STOP.\n');
        end 
        %% Update
        g = g+1;
    end
    if PRINT_SHOW==1
        figure; hold on; title('(mu/mu_{I},lambda)-sigmaSA-ES');
            plot(f_g, 'DisplayName','f');
            plot(sigma_g, 'DisplayName','sigma');       
            % plot(r_g, 'DisplayName','R');    
            yscale('log'); legend;
    end

end