function out = prob_normal_distribution(a,b_sq)

    out = (1/sqrt(2*pi*b_sq))*exp((-0.5)*a^2/b_sq);
    
end