function out = prob_triangular_distribution(a,b_sq)

    out = max(0,1/sqrt(6*b_sq)-abs(a)/(5*b_sq));
    
end