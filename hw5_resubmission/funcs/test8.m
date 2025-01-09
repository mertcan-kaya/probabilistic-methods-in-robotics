clc, clear all, close all

app.nSmp = 10000;

app.prob_type = 0;
app.traj_type = 1;

app.dt = 1;
app.dist_type = 0;

if app.prob_type == 0
    app.alpha1 = 0.001;
    app.alpha2 = 0.001;
    app.alpha3 = 0.001;
    app.alpha4 = 0.001;
    app.alpha5 = 0.001;
    app.alpha6 = 0.001;
%     app.alpha1 = 0.0;
%     app.alpha2 = 0.0;
%     app.alpha3 = 0.0;
%     app.alpha4 = 0.0;
%     app.alpha5 = 0.0;
%     app.alpha6 = 0.0;
else
    app.alpha1 = 0.00001;
    app.alpha2 = 0.00001;
    app.alpha3 = 0.00001;
    app.alpha4 = 0.00001;
end

sigma_d = 0.025;
sigma_phi = deg2rad(2.5);
sigma_s = 0.4;

%% generate trajectory
if app.traj_type == 0
    nStp = 8;
    dt = app.dt;

    x = 1;
    y = 0;
    theta = -pi/2;

    v = pi/4*ones(1,nStp); % m/s
    omega = -pi/4*ones(1,nStp); % rad/s

    x_array     = x*ones(1,nStp);
    y_array     = y*ones(1,nStp);
    theta_array = theta*ones(1,nStp);

    for i = 1:nStp-1
        x_array(i+1) = x_array(i) - (v(i)/omega(i))*(sin(theta_array(i)) - sin(theta_array(i)+omega(i)*app.dt));
        y_array(i+1) = y_array(i) + (v(i)/omega(i))*(cos(theta_array(i)) - cos(theta_array(i)+omega(i)*app.dt));
        theta_array(i+1) = normalizeFunc(theta_array(i) + omega(i)*app.dt,-pi,pi);
    end

else
    dt = 1;
    
    nStp = 3;

    x_array = [0,2,4,6,8,8,8,8,8,8,8,6,4,2,0];
    y_array = [0,0,0,0,0,0,2,4,6,8,8,8,8,8,8];
    theta_array = [0,0,0,0,0,pi/2,pi/2,pi/2,pi/2,pi/2,pi,pi,pi,pi,pi];
        
    v = zeros(1,nStp); % m/s
    omega = zeros(1,nStp); % rad/s

    for i = 1:nStp-1
        v(1,i) = sqrt((x_array(i+1)-x_array(i))^2+(y_array(i+1)-y_array(i))^2); % m/s
        omega(1,i) = theta_array(i+1)-theta_array(i); % rad/s
    end

end
x_t_array = [x_array;y_array;theta_array];
v_t_array = [v;omega];

% landmark
point(1).x = x_array(1,1);
point(1).y = y_array(1,1);
point(1).s = 'o';

point(2).x = x_array(1,2);
point(2).y = y_array(1,2);
point(2).s = '^';

point(3).x = x_array(1,3);
point(3).y = y_array(1,3);
point(3).s = '^';

nP = length(point);

xM = zeros(1,nP);
yM = zeros(1,nP);
sM = zeros(1,nP);
for i = 1:nP
    xM(i) = point(i).x;
    yM(i) = point(i).y;
    sM(i) = point(i).s;
end

app.L_info = [xM;yM;sM];


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

% figure
% plot(x_t_array(1,:),x_t_array(2,:),'.k')
% ylim([min(x_t_array(2,:))-0.25 max(x_t_array(2,:))+0.25])
% hold('on')
% cmap = jet(nStp);
% x_steps = zeros(3,app.nSmp,nStp);
% h = zeros(nStp,1);
% for i = 1:nStp
%     x_steps(:,:,i) = x_t(:,i:nStp:app.nSmp*nStp);
%     h(i) = plot(x_t(1,i:nStp:app.nSmp*nStp),x_t(2,i:nStp:app.nSmp*nStp),'.','Color', cmap(i, :));
% end
% axis equal
% grid on

%% EKF localization known correspondences

mu_bar = zeros(3,1,nStp);
mu_bar(:,:,1) = x_t_array(:,1);
Sigma_bar = zeros(3,3,nStp);

for k = 1:nStp-1
    if app.prob_type == 0
        u_t = v_t_array(:,k);
    else
        u_t = [x_t_array(:,k);x_t_array(:,k+1)];
    end
    mu_bar(:,:,k+1) = g_t(app,u_t,mu_bar(:,:,k));
    Sigma_bar(:,:,k+1) = G_t(app,u_t,mu_bar(:,:,k))*Sigma_bar(:,:,k)*G_t(app,u_t,mu_bar(:,:,k))'+V_t(app,u_t,mu_bar(:,:,k))*M_t(app,u_t)*V_t(app,u_t,mu_bar(:,:,k))';

    % Correction
    Q_t = diag([sigma_d,sigma_phi,sigma_s]);
    mu_t_bar = mu_bar(:,:,k+1);
    Sigma_t_bar = Sigma_bar(:,:,k+1);
    for i = 1:nP
        q = (m_j(1)-mu_t_bar(1))^2+(m_j(2)-mu_t_bar(2))^2;
        z_t_hat = [sqrt(q);atan2(m_j(2)-mu_t_bar(2),m_j(1)-mu_t_bar(1))-mu_t_bar(3);m_j(3)];
        H_t = [ -(m_j(1)-mu_t_bar(1))/sqrt(q), -(m_j(2)-mu_t_bar(2))/sqrt(q), 0
                (m_j(2)-mu_t_bar(2))/q, -(m_j(1)-mu_t_bar(1))/q, -1
                0, 0, 0];
        S_t = H_t*Sigma_t_bar*H_t'+Q_t;
        K_t = Sigma_t_bar*H_t'/S_t;
        z_t = h_t(x_t_array(:,k),m_j);
        mu_t_bar = mu_t_bar + K_t*(z_t-z_t_hat);
        Sigma_t_bar = (eye(3)-K_t*H_t)*Sigma_t_bar;
    end
end

% figure
% grid on
% hold on
% axis equal
% mu = zeros(3,nStp);
% Sigma = zeros(3,3,nStp);
% for k = 1:nStp
%     mu(:,k) = mean(x_steps(:,:,k),2);
%     S12 = ((x_steps(1,:,k)-mu(1,k))*(x_steps(2,:,k)-mu(2,k))')/(app.nSmp-1);
%     S13 = ((x_steps(1,:,k)-mu(1,k))*(x_steps(3,:,k)-mu(3,k))')/(app.nSmp-1);
%     S23 = ((x_steps(2,:,k)-mu(2,k))*(x_steps(3,:,k)-mu(3,k))')/(app.nSmp-1);
%     Sigma(:,:,k) = [var(x_steps(1,:,k)) S12 S13;S12 var(x_steps(2,:,k)) S23;S13 S23 var(x_steps(3,:,k))];
%     
%     drawEllipse(mu(:,k),Sigma(:,:,k),cmap(k, :))
% end
% 
% R = Sigma(:,:,2)

%% Extended Kalman Filter

mu_bar_e = zeros(3,nStp-1);
Sigma_bar_e = zeros(3,3,nStp-1);

% prediction
mu_bar_e(:,1)       = x_t_array(:,1);
Sigma_bar_e(:,:,1)	= Sigma(:,:,1);

drawEllipse(mu_bar_e(:,1),Sigma_bar_e(:,:,1),'#000000')

for k = 1:nStp-1
    if app.prob_type == 0
        u_t = v_t_array(:,k);
    else
      	u_t = [x_t_array(:,k);x_t_array(:,k+1)];
    end
    g_tk = g_t(app,u_t,mu_bar_e(:,k));
    G_tk = G_t(app,u_t,mu_bar_e(:,k));

    % prediction
    mu_bar_e(:,k+1)       = g_tk;
    Sigma_bar_e(:,:,k+1)	= G_tk*Sigma_bar_e(:,:,k)*G_tk' + R;

    drawEllipse(mu_bar_e(:,k+1),Sigma_bar_e(:,:,k+1),'#000000')
end
plot(x_t_array(1,:),x_t_array(2,:),'.k')

%% Unscented Kalman Filter

n = 3;
chiNb = 2*n + 1;

kappa = 10;      % kappa >= 0
alpha = 0.5;    % alpha € (0,1]
lambda = alpha^2*(n+kappa)-n;
beta = 2;

mu_bar_u = zeros(3,1,nStp);
Sigma_bar_u = zeros(3,3,nStp);

mu_bar_u(:,:,1) = x_t_array(:,1);
Sigma_bar_u(:,:,1) = Sigma(:,:,1);

drawEllipse(mu_bar_u(:,:,1),Sigma_bar_u(:,:,1),'#FF00FF')
    
chi_bar = zeros(3,chiNb,nStp-1);
for k = 1:nStp-1
    chi = zeros(3,chiNb);
    mu_offset = sqrtm((n+lambda)*Sigma_bar_u(:,:,k));
    for j = 1:chiNb
        if j == 1
            chi(:,j) = mu_bar_u(:,:,k);
        elseif j <= n+1
            chi(:,j) = mu_bar_u(:,:,k) + mu_offset(:,j-1);
        else
            chi(:,j) = mu_bar_u(:,:,k) - mu_offset(:,j-n-1);
        end
    end

    omega_m = zeros(1,chiNb);
    omega_c = zeros(1,chiNb);
    for j = 1:chiNb
        if j == 1
            omega_m(j) = lambda/(n+lambda);
            omega_c(j) = lambda/(n+lambda)+1-alpha^2+beta;
        else
            omega_m(j) = 1/(2*(n+lambda));
            omega_c(j) = 1/(2*(n+lambda));
        end
    end

    chi_tk_bar_s = zeros(3,chiNb);
    if app.prob_type == 0
        u_t = v_t_array(:,k);
    else
      	u_t = [x_t_array(:,k);x_t_array(:,k+1)];
    end
    for j = 1:chiNb
        chi_tk_bar_s(:,j) = g_t(app,u_t,chi(:,j));
    end

    for j = 1:chiNb
        mu_bar_u(:,:,k+1) = mu_bar_u(:,:,k+1) + omega_m(j)*chi_tk_bar_s(:,j);
    end

    for j = 1:chiNb
        Sigma_bar_u(:,:,k+1) = Sigma_bar_u(:,:,k+1) + omega_c(j)*(chi_tk_bar_s(:,j)-mu_bar_u(:,:,k+1))*(chi_tk_bar_s(:,j)-mu_bar_u(:,:,k+1))';
    end
    Sigma_bar_u(:,:,k+1) = Sigma_bar_u(:,:,k+1) + R;

    mu_tk_offset_bar = sqrtm((n+lambda)*Sigma_bar_u(:,:,k+1));
    for j = 1:chiNb
        if j == 1
            chi_bar(:,j,k+1) = mu_bar_u(:,:,k+1);
        elseif j <= n+1
            chi_bar(:,j,k+1) = mu_bar_u(:,:,k+1) + mu_tk_offset_bar(:,j-1);
        else
            chi_bar(:,j,k+1) = mu_bar_u(:,:,k+1) - mu_tk_offset_bar(:,j-n-1);
        end
    end
    
    drawEllipse(mu_bar_u(:,:,k+1),Sigma_bar_u(:,:,k+1),'#FF00FF')
    plot(chi_bar(1,[1:3 5:6],k+1),chi_bar(2,[1:3 5:6],k+1),'Marker','.','LineStyle','none','Color','#FF00FF','LineWidth',1);

end

% for k = 2:nStp
%     figure
%     hold on
%     % plot(x_t_array(1,:),x_t_array(2,:),'.k')
%     drawEllipse(mu(:,k),Sigma(:,:,k),cmap(k, :))
%     drawEllipse(mu_bar(:,k),Sigma_bar(:,:,k),'#000000')
%     drawEllipse(mu_bar_u(:,:,k),Sigma_bar_u(:,:,k),'#FF00FF')
%     plot(chi_bar(1,[1:3 5:6],k),chi_bar(2,[1:3 5:6],k),'Marker','.','LineStyle','none','Color','#FF00FF','LineWidth',1);
% end

figure(5)
grid on
hold on
axis equal
for k = 2:nStp-1
    drawEllipse(mu(:,k),Sigma(:,:,k),cmap(k, :))
    drawEllipse(mu_bar_e(:,k),Sigma_bar_e(:,:,k),'c')
    drawEllipse(mu_bar(:,:,k),Sigma_bar(:,:,k),'k')
    drawEllipse(mu_t_bar,Sigma_t_bar,'m')
end

function drawEllipse(mu,Sigma,color)

    if sum(isnan(Sigma),'all')
        return
    else
%         the 95% confidence interval error ellipse
        chisquare_val = sqrt(5.991);
        
%         % the 99% confidence interval error ellipse
%         chisquare_val = sqrt(9.210);
        
        theta_grid = linspace(0,2*pi);

        [eigenvec, eigenval] = eig(Sigma(1:2,1:2));

        % Get the index of the largest eigenvector
        [largest_eigenvec_ind_c, ~] = find(eigenval == max(max(eigenval)));
        largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

        % Get the largest eigenvalue
        largest_eigenval = max(max(eigenval));

        % Get the smallest eigenvector and eigenvalue
        if(largest_eigenvec_ind_c == 1)
            smallest_eigenval = max(eigenval(:,2));
%             smallest_eigenvec = eigenvec(:,2);
        else
            smallest_eigenval = max(eigenval(:,1));
%             smallest_eigenvec = eigenvec(1,:);
        end

        % Calculate the angle between the x-axis and the largest eigenvector
        angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

        % This angle is between -pi and pi.
        % Let's shift it such that the angle is between 0 and 2pi
        if(angle < 0)
            angle = angle + 2*pi;
        end

        phi = angle;
        a = chisquare_val * sqrt(largest_eigenval);
        b = chisquare_val * sqrt(smallest_eigenval);

        % the ellipse in x and y coordinates 
        ellipse_x_r  = a*cos(theta_grid);
        ellipse_y_r  = b*sin(theta_grid);

        % rotation matrix of angle phi
        rot = [cos(phi) sin(phi); -sin(phi) cos(phi)];

        % rotation of the ellipse to angle phi
        r_ellipse = [ellipse_x_r;ellipse_y_r]' * rot;

        plot(mu(1),mu(2),'Marker','*','Color',color);
        plot(r_ellipse(:,1) + mu(1),r_ellipse(:,2) + mu(2),'Color',color);
    end
    
end

function n = normalizeFunc(value,start,finish) 

  width       = finish - start   ;
  offsetValue = value - start ;

  n = ( offsetValue - ( floor( offsetValue / width ) * width ) ) + start ;
  
end

function x_t = g_t(app,u_t,x_tprev)

    theta   = x_tprev(3);
    
    if app.prob_type == 0
        v       = u_t(1);
        omega   = u_t(2);
    
        if abs(omega) > 1e-5
            x_t = x_tprev + [  -(v/omega)*(sin(theta) - sin(theta+omega*app.dt))
                                (v/omega)*(cos(theta) - cos(theta+omega*app.dt))
                                omega*app.dt ];
        else
            x_t = x_tprev + [	v*cos(theta + omega*app.dt)
                              	v*sin(theta + omega*app.dt)
                             	omega*app.dt ];
        end
    else
        x_bar       = u_t(1);
        y_bar       = u_t(2);
        theta_bar	= u_t(3);

        x_bar_prm       = u_t(4);
        y_bar_prm       = u_t(5);
        theta_bar_prm   = u_t(6);
    
        x_t = x_tprev + [	sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*cos(theta + atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar)
                            sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*sin(theta + atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar)
                            theta_bar_prm - theta_bar ];
    end

end

function g_t_prm = G_t(app,u_t,x_tprev)

    theta   = x_tprev(3);
    
    if app.prob_type == 0
        v       = u_t(1);
        omega   = u_t(2);
    
        if abs(omega) > 1e-5
            g_t_prm = [	1,0,(v/omega)*(cos(omega*app.dt + theta) - cos(theta))
                        0,1,(v/omega)*(sin(omega*app.dt + theta) - sin(theta))
                        0,0,1 ];
        else
            g_t_prm = [	1,0,-v*sin(omega*app.dt + theta)
                        0,1,v*cos(omega*app.dt + theta)
                        0,0,1 ];
        end
    else
        x_bar       = u_t(1);
        y_bar       = u_t(2);
        theta_bar	= u_t(3);

        x_bar_prm 	= u_t(4);
        y_bar_prm 	= u_t(5);
    
        g_t_prm = [	1, 0, -sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*sin(atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar + theta)
                  	0, 1,  sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*cos(atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar + theta)
                   	0, 0,  1 ];
    end

end

function out = V_t(app,u_t,x_tprev)

    theta   = x_tprev(3);
    
    v       = u_t(1);
    omega   = u_t(2);

    if abs(omega) > 1e-5
        out = [	(-sin(theta)+sin(theta+omega*app.dt))/omega, v*((sin(theta)-sin(theta+omega*app.dt))/omega+cos(theta+omega*app.dt)*app.dt)/omega, 0
                (cos(theta)-cos(theta+omega*app.dt))/omega, v*(-(cos(theta)-cos(theta+omega*app.dt))/omega+sin(theta+omega*app.dt)*app.dt)/omega, 0
                0                                        , app.dt                                                                               , app.dt];
    else

        out = [ cos(theta+omega*app.dt), -v*sin(theta+omega*app.dt)*app.dt, 0
                sin(theta+omega*app.dt),  v*cos(theta+omega*app.dt)*app.dt, 0
                0                       ,  app.dt                           , app.dt];
    end
        
end

function out = M_t(app,u_t)

    v       = u_t(1);
    omega   = u_t(2);

    out = [ app.alpha1*v^2+app.alpha2*omega^2, 0                                , 0
            0                               ,  app.alpha3*v^2+app.alpha4*omega^2, 0
            0                               ,  0                                , app.alpha5*v^2+app.alpha6*omega^2];
        
end

function out = h_t(x_t,m_j)

    out = [ sqrt((m_j(1)-x_t(1))^2+(m_j(2)-x_t(2))^2)
            atan2(m_j(2)-x_t(2),m_j(1)-x_t(1))-x_t(3)
            m_j(3) ];
        
end
