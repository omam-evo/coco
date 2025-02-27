function [fmed_g] = stag_check(fmed_g,arfitness,mu,g)
    if g>10
        fmed_g(1:end-1) = fmed_g(2:end);
        fmed_g(end) = median(arfitness(1:mu));
    else
        fmed_g(g) = median(arfitness(1:mu));
    end
end

