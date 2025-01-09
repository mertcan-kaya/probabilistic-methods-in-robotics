function out = sample_triangular_distribution(b_sq)

    b = sqrt(b_sq);
    
    out = (sqrt(6)/2)*(rand2(-b,b) + rand2(-b,b));
    
end