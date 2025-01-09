function x_t = hypothesized_motion_model_velocity(app,u_t,x_tprev)

    v       = u_t(1);
    omega   = u_t(2);
    
    x       = x_tprev(1);
    y       = x_tprev(2);
    theta   = x_tprev(3);
    
    if abs(omega) > 1e-5
        x_prm       = x - (v/omega)*(sin(theta) - sin(theta + omega*app.dt));
        y_prm       = y + (v/omega)*(cos(theta) - cos(theta + omega*app.dt));
    else
        x_prm       = x + v*cos(theta + omega*app.dt);
        y_prm       = y + v*sin(theta + omega*app.dt);
    end
    theta_prm   = theta + omega*app.dt;
    
    x_t = [x_prm;y_prm;theta_prm];
    
end