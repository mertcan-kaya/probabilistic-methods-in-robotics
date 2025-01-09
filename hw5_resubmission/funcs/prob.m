function out = prob(app,a,b_sq)

    if app.dist_type == 1
        out = prob_triangular_distribution(a,b_sq);
    else
        out = prob_normal_distribution(a,b_sq);
    end
    
end