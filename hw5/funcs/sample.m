function out = sample(app,b_sq)

    if app.dist_type == 1
        out = sample_triangular_distribution(b_sq);
    else
        out = sample_normal_distribution(b_sq);
    end
    
end