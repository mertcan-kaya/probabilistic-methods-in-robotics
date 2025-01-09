function x_t = sample_motion_model_velocity(app,u_t,x_tprev)

    v       = u_t(1);
    omega   = u_t(2);
    
    x       = x_tprev(1);
    y       = x_tprev(2);
    theta   = x_tprev(3);
    
    b_sq1   = app.alpha1*v^2 + app.alpha2*omega^2;
    b_sq2   = app.alpha3*v^2 + app.alpha4*omega^2;
    b_sq3   = app.alpha5*v^2 + app.alpha6*omega^2;
    
    v_hat       = v + sample(app,b_sq1);
    omega_hat   = omega + sample(app,b_sq2);
    gamma_hat   = sample(app,b_sq3);
    
    if abs(omega_hat) > 1e-5
        x_prm       = x - (v_hat/omega_hat)*(sin(theta) - sin(theta + omega_hat*app.dt));
        y_prm       = y + (v_hat/omega_hat)*(cos(theta) - cos(theta + omega_hat*app.dt));
        theta_prm   = theta + (omega_hat + gamma_hat)*app.dt;
    else
        x_prm       = x + v_hat*cos(theta + omega_hat*app.dt);
        y_prm       = y + v_hat*sin(theta + omega_hat*app.dt);
        theta_prm   = theta + (omega_hat + gamma_hat)*app.dt;
    end
    
    x_t = [x_prm;y_prm;theta_prm];
    
end