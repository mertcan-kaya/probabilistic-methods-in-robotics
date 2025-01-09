function out = sample_normal_distribution(b_sq)

    b = sqrt(b_sq);
    
    sumrand = 0;
    for i = 1:12
        sumrand = sumrand + rand2(-b,b);
    end
    
    out = 0.5*sumrand;
    
end