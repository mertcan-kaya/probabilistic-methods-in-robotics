function x_t = sample_motion_model_odometry(app,u_t,x_tprev)

    x_bar       = u_t(1);
    y_bar       = u_t(2);
    theta_bar	= u_t(3);
    
    x_bar_prm       = u_t(4);
    y_bar_prm       = u_t(5);
    theta_bar_prm   = u_t(6);
    
    x       = x_tprev(1);
    y       = x_tprev(2);
    theta   = x_tprev(3);
    
    delta_rot1  = atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar;
    delta_trans = sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2);
    delta_rot2  = theta_bar_prm - theta_bar - delta_rot1;
    
    b_sq1 = app.alpha1*delta_rot1^2 + app.alpha2*delta_trans^2;
    b_sq2 = app.alpha3*delta_trans^2 + app.alpha4*(delta_rot1^2+delta_rot2^2);
    b_sq3 = app.alpha1*delta_rot2^2 + app.alpha2*delta_trans^2;
    
    delta_rot1_hat	= delta_rot1 - sample(app,b_sq1);
    delta_trans_hat = delta_trans - sample(app,b_sq2);
    delta_rot2_hat	= delta_rot2 - sample(app,b_sq3);
    
    x_prm       = x + delta_trans_hat*cos(theta + delta_rot1_hat);
    y_prm       = y + delta_trans_hat*sin(theta + delta_rot1_hat);
    theta_prm   = theta + delta_rot1_hat + delta_rot2_hat;
    
    x_t = [x_prm;y_prm;theta_prm];
    
end