clc, clear all, close all

app.nSmp = 100;

app.prob_type = 0;
app.traj_type = 1;
app.comp_type = 1;

app.dt = 1;
app.dist_type = 0;

if app.prob_type == 0
    app.alpha1 = 0.001;
    app.alpha2 = 0.001;
    app.alpha3 = 0.001;
    app.alpha4 = 0.001;
    app.alpha5 = 0.001;
    app.alpha6 = 0.001;
else
    app.alpha1 = 0.01;
    app.alpha2 = 0.01;
    app.alpha3 = 0.01;
    app.alpha4 = 0.01;
end

%% generate trajectory
if app.traj_type == 0
    nStp = 8;

    x = 1;
    y = 0;
    theta = -pi/2;

    v = pi/4*ones(1,nStp); % m/s
    omega = -pi/4*ones(1,nStp); % rad/s

    x_array     = x*ones(1,nStp);
    y_array     = y*ones(1,nStp);
    theta_array = theta*ones(1,nStp);

    for i = 1:nStp-1
        x_array(i+1) = x_array(i) - (v(i)/omega(i))*(sin(theta_array(i)) - sin(theta_array(i)+omega(i)));
        y_array(i+1) = y_array(i) + (v(i)/omega(i))*(cos(theta_array(i)) - cos(theta_array(i)+omega(i)));
        theta_array(i+1) = theta_array(i) + omega(i);
    end

else
    if app.prob_type == 0
        nStp = 15;

        x_array = [0,2,4,6,8,8,8,8,8,8,8,6,4,2,0];
        y_array = [0,0,0,0,0,0,2,4,6,8,8,8,8,8,8];
        theta_array = [0,0,0,0,0,pi/2,pi/2,pi/2,pi/2,pi/2,pi,pi,pi,pi,pi];
    else
        nStp = 13;

        x_array = [0,2,4,6,8,8,8,8,8,6,4,2,0];
        y_array = [0,0,0,0,0,2,4,6,8,8,8,8,8];
        theta_array = [0,0,0,0,0,pi/2,pi/2,pi/2,pi,pi,pi,pi,pi];
    end

    v = zeros(1,nStp); % m/s
    omega = zeros(1,nStp); % rad/s

    for i = 1:nStp-1
        v(1,i) = sqrt((x_array(i+1)-x_array(i))^2+(y_array(i+1)-y_array(i))^2); % m/s
        omega(1,i) = theta_array(i+1)-theta_array(i); % rad/s
    end

end
x_t_array = [x_array;y_array;theta_array];
v_t_array = [v;omega];


if app.comp_type == 0
    % sampling
    x_t = zeros(3,app.nSmp*nStp);
    for j = 1:app.nSmp

        k = 1;
        x_t(:,1+(j-1)*nStp) = x_t_array(:,k);

        for i = 1+(j-1)*nStp:j*(nStp)-1

            if app.prob_type == 0
                u_t = v_t_array(:,k);
                x_t(:,i+1) = sample_motion_model_velocity(app,u_t,x_t(:,i));
            else
                u_t = [x_t_array(:,k);x_t_array(:,k+1)];
                x_t(:,i+1) = sample_motion_model_odometry(app,u_t,x_t(:,i));
            end

            k = k+1;
        end
    end

    figure
    plot(x_t_array(1,:),x_t_array(2,:),'.k')
    ylim([min(x_t_array(2,:))-0.25 max(x_t_array(2,:))+0.25])
    hold('on')
    cmap = jet(nStp);
    h = zeros(nStp,1);
    for i = 1:nStp
        h(i) = plot(x_t(1,i:nStp:app.nSmp*nStp),x_t(2,i:nStp:app.nSmp*nStp),'.','Color', cmap(i, :));
        if i == 1
            set(h(i),'Marker','x')
        end
    end
    % ylim('auto')
    ylim([min(x_t(1,:))-0.25 max(x_t(2,:))+0.25])
    axis equal

else
    % direct evaluation

    if app.traj_type == 0
        x_grid = -2:0.01:2;
        y_grid = -2:0.01:2;
    else
        x_grid = -1:0.01:10;
        y_grid = -1:0.01:10;
    end

    p_x_t = zeros(length(x_grid),length(y_grid),nStp);
    for k = 1:nStp-1

        for j = 1:length(x_grid)
            for i = 1:length(y_grid)
                x_t = [x_grid(j);y_grid(i);0];
                if app.prob_type == 0
                    u_t = v_t_array(:,k);
                    p_x_t(i,j,k) = motion_model_velocity(app,x_t,u_t,x_t_array(:,k));
                else
                    u_t = [x_t_array(:,k);x_t_array(:,k+1)];
                    p_x_t(i,j,k) = motion_model_odometry(app,x_t,u_t,x_t_array(:,k));
                end
            end
        end

    end

    figure
    plot(x_t_array(1,:),x_t_array(2,:),'.k')
    ylim([min(x_t_array(2,:))-0.25 max(x_t_array(2,:))+0.25])
    axis equal
    % figure
    % surf(x_grid,y_grid,p_x_t_sum)
    % shading interp

    dim = 3;

    for k = 1:nStp-1
        figure
        plot(x_t_array(1,:),x_t_array(2,:),'.k')
        hold on
        azSurf1 = surf(x_grid,y_grid,p_x_t(:,:,k));
        xlabel('x')
        ylabel('y')
        shading interp
        set(azSurf1, 'FaceAlpha', 'interp','AlphaData',p_x_t(:,:,k));
        % axis equal
        if app.traj_type == 0
            zlim([0 9])
        else
            zlim([0 5])
        end
        view(dim)
    end
    
end