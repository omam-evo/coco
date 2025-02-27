function [res] = terminate_cma(bbob,STOP,fmin,g,sigma,max_D,min_D,x_new,x_old,f_med_g)
    res = 0;
    if ~isempty(bbob)         % if script is called with bbob-problem
        if cocoProblemFinalTargetHit(bbob) %cocoProblemGetBestValue(problem)
            fprintf('SUCCESS cma.\n');
            res = 1;
        end
    elseif fmin <= STOP.F_STOP % if script is called with F_STOP (no bbob)
        fprintf('F STOP cma.\n');
        res = 1;
    end
    if g>=STOP.G_MAX
        fprintf('GEN STOP cma.\n');
        res = 1;
    elseif sigma*max_D < STOP.SIGMA_D_STOP % sigma and maximum eigenvalue stop
        fprintf('SIGMA*D STOP cma.\n');
        res = 1;
    elseif max_D/min_D > STOP.COND_SQRT   % sqrt of condition 1e14 => 1e7
        fprintf('COND STOP cma.\n');
        res = 1;
    elseif norm(x_new-x_old) < STOP.X_TOL % update
        fprintf('XTOL STOP cma.\n');
        res = 1;  
    elseif g>10 && max(f_med_g)-min(f_med_g) < STOP.F_STAG_TOL % stagnation
        fprintf('F_STAG_TOL STOP cma.\n');
        res = 1; 
    end
end

