function prob_out = motion_model_odometry(app,x_t,u_t,x_tprev)

    x_prm       = x_t(1);
    y_prm       = x_t(2);
    theta_prm	= x_t(3);
    
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
    
    delta_rot1_hat	= atan2(y_prm-y, x_prm-x) - theta;
    delta_trans_hat = sqrt((x-x_prm)^2 + (y-y_prm)^2);
    delta_rot2_hat	= theta_prm - theta - delta_rot1_hat;
    
    err_rot1 = delta_rot1-delta_rot1_hat;
    err_trans = delta_trans-delta_trans_hat;
    err_rot2 = delta_rot2-delta_rot2_hat;
    
    b_sq1 = app.alpha1*delta_rot1_hat^2 + app.alpha2*delta_trans_hat^2;
    b_sq2 = app.alpha3*delta_trans_hat^2 + app.alpha4*(delta_rot1_hat^2+delta_rot2_hat^2);
    b_sq3 = app.alpha1*delta_rot2_hat^2 + app.alpha2*delta_trans_hat^2;
    
%     if b_sq1 > 1e-7
        prob_rot1 = prob(app,err_rot1,b_sq1);
%     else
%         prob_rot1 = 1;
%     end
%     if b_sq2 > 1e-7
        prob_trans = prob(app,err_trans,b_sq2);
%     else
%         prob_trans = 1;
%     end
%     if b_sq3 > 1e-7
%         prob_rot2 = prob(app,err_rot2,b_sq3);
%     else
%         prob_rot2 = 1;
%     end
    
    prob_out = prob_rot1*prob_trans;
    
end