function prob_out = motion_model_velocity(app,x_t,u_t,x_tprev)

    x_prm       = x_t(1);
    y_prm       = x_t(2);
    theta_prm	= x_t(3);
    
    v       = u_t(1);
    omega   = u_t(2);
    
    x       = x_tprev(1);
    y       = x_tprev(2);
    theta   = x_tprev(3);
    
    mu = 0.5*( ( (x-x_prm)*cos(theta)+(y-y_prm)*sin(theta) ) / ( (y-y_prm)*cos(theta)-(x-x_prm)*sin(theta) ) );
    
    x_str = (x+x_prm)/2 + mu*(y-y_prm);
    y_str = (y+y_prm)/2 + mu*(x_prm-x);
    r_str = sqrt( (x-x_str)^2 + (y-y_str)^2 );
%     r_str2 = sqrt( (x_prm-x_str)^2 + (y_prm-y_str)^2 )
    
    atan2(y_prm-y_str, x_prm-x_str);
    atan2(y-y_str, x-x_str);
    dtheta = atan2(y_prm-y_str, x_prm-x_str) - atan2(y-y_str, x-x_str);
    
    v_hat       = -(dtheta/app.dt)*r_str;
    omega_hat   = dtheta/app.dt;
    gamma_hat   = ( (theta_prm-theta)/app.dt ) - omega_hat;
    
    err_v = v-v_hat;
    err_omega = omega-omega_hat;
    err_gamma = gamma_hat;
        
    b_sq1   = app.alpha1*v^2 + app.alpha2*omega^2;
    b_sq2   = app.alpha3*v^2 + app.alpha4*omega^2;
    b_sq3   = app.alpha5*v^2 + app.alpha6*omega^2;
    
%     if b_sq1 > 1e-7
        prob_v = prob(app,err_v,b_sq1);
%     else
%         prob_v = 1;
%     end
%     if b_sq2 > 1e-7
        prob_omega = prob(app,err_omega,b_sq2);
%     else
%         prob_omega = 1;
%     end
%     if b_sq3 > 1e-7
        prob_gamma = prob(app,err_gamma,b_sq3);
%     else
%         prob_gamma = 1;
%     end
    
    prob_out = prob_v*prob_omega*prob_gamma;
    
%     a = ((v-v_hat)*b_sq2 + (omega-omega_hat)*b_sq1)/(b_sq2+b_sq1);
%     b_sq = (b_sq1*b_sq2)/(b_sq1+b_sq2);
    
end