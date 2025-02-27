function [xmin,fmin,counteval] = bipop_cma_es(fun, bbob, mu0, lambda0, x0, x_low, x_up, sigma0, STOP, PRINT_SHOW)

    PRINT_SHOW_INNER = 0;
    N = length(x0);
    STOP_UPD = STOP; % needed to update current number of feval
    
    [xmin,fmin,counteval,~,~,~,~,~] = ...
        purecmaes_basic(fun, bbob, mu0, lambda0, x0, sigma0, STOP, PRINT_SHOW_INNER);

    % Terminate if simple run was successful
    if ~isempty(bbob)
        if cocoProblemFinalTargetHit(bbob) % cocoProblemGetBestValue(problem)
            fprintf('SUCCESS bipop.\n');
            return;
        end
    elseif fmin <= STOP.F_STOP % not needed for coco
        fprintf('F STOP bipop.\n');
        return;
    end

    b1 = counteval;
    n = 0;
    n_small = 0;
    b_large = 0;
    b_small = 0;
    lambda_max = 2^9*lambda0; % z.B. N=40: LAM=7680
    lambda = lambda0;

    while counteval < STOP.FEVAL_MAX && lambda*2<=lambda_max
        
        n = n+1;
        lambda =  2^n*lambda0; %HGB: 2^(n-n_small)*lambda0
        x = rand(N,1).*(x_up-x_low)+x_low;
        STOP_UPD.FEVAL_MAX = STOP.FEVAL_MAX-counteval; % remaining number of feval

        if n>1 && b_small<b_large % Small pop. regime
  
            %STOP_UPD.FEVAL_MAX = STOP_UPD.FEVAL_MAX/2;  %Hansen: use only half of the remaining budget for local search
            sigma_0_small = sigma0/(100^rand);
            lambda_small = floor(lambda0*(0.5*lambda/lambda0)^(rand^2));
            
            [xres,fres,feval_small,~,~,~,~,~] = ...
                purecmaes_basic(fun, bbob, floor(lambda_small/2), lambda_small, x, sigma_0_small, STOP_UPD, PRINT_SHOW_INNER);
            b_small = b_small + feval_small;
            n_small = n_small + 1;
            fprintf(['\t .. small:', ' lambda:' num2str(lambda_small,'%i'), ' b_small:' num2str(b_small), ' fres:' num2str(fres),'\n']);
        else % Large pop. regime
            
            [xres,fres,feval_large,~,~,~,~,~] = ...
                purecmaes_basic(fun, bbob, floor(lambda/2), lambda, x, sigma0, STOP_UPD, PRINT_SHOW_INNER);   
            b_large = b_large + feval_large;
            fprintf(['\t .. large:', ' lambda:' num2str(lambda,'%i'), ' b_large:' num2str(b_large), ' fres:' num2str(fres),'\n']);
        end
        
        counteval = b1 + b_small + b_large;
        
        if fres<fmin
            xmin = xres;
            fmin = fres;
        end

        %% Output 
        if PRINT_SHOW==1
            disp(['n:',num2str(n,'%i'), ' feval:', num2str(counteval,'%i'), ' fmin:' num2str(fmin)]);
            f_g(n) = fmin;
        end

        %% Terminate outer loop if target has been reached
        % continue loop if func. eval. budget remaining
        if ~isempty(bbob)
            if cocoProblemFinalTargetHit(bbob) %cocoProblemGetBestValue(problem)
                fprintf('SUCCESS bipop.\n');
                break
            end
        elseif fmin <= STOP.F_STOP % not needed for coco
            fprintf('F STOP bipop.\n');
            break;
        end

    end % function evaluations
    if counteval >= STOP.FEVAL_MAX
        fprintf('FEVAL STOP bipop.\n');
    end
    if PRINT_SHOW==1
        figure; hold on; title('bipop-cma-es');
            plot(f_g, 'DisplayName','f');
            yscale('log'); legend;
    end

end % function
