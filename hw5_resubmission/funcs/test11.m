clc, clear all, close all

app.nSmp = 10000;

app.prob_type = 0;
app.traj_type = 0;

app.dt = 1;
app.dist_type = 0;

if app.prob_type == 0
    app.alpha1 = 0.01;
    app.alpha2 = 0.01;
    app.alpha3 = 0.01;
    app.alpha4 = 0.01;
    app.alpha5 = 0.01;
    app.alpha6 = 0.01;
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

app.sigma_d = 0.005;
app.sigma_phi = deg2rad(0.5);
app.sigma_s = 0.01;

% app.sigma_d = 0.000;
% app.sigma_phi = deg2rad(0);
% app.sigma_s = 0.4;

%% generate trajectory
if app.traj_type == 0
    nStp = 4;
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
% point(1).x = x_array(1,1)+0.2;
% point(1).y = y_array(1,1)-0.4;
% point(1).s = 'o';

point(1).x = x_array(1,2)+0.1;
point(1).y = y_array(1,2)-0.45;
point(1).s = 'o';

point(2).x = x_array(1,3)-0.2;
point(2).y = y_array(1,3)+0.2;
point(2).s = '^';

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

%% EKF localization known correspondences

mu_bar = zeros(3,1,nStp);
Sigma_bar = zeros(3,3,nStp);

mu = zeros(3,1,nStp);
Sigma = zeros(3,3,nStp);

mu_bar(:,:,1) = x_t_array(:,1);
Sigma_bar(:,:,1) = zeros(3,3);

mu(:,:,1) = mu_bar(:,:,1);
Sigma(:,:,1) = Sigma_bar(:,:,1);

m = zeros(3,nP);
m_hat = zeros(3,nP);
z_t = zeros(3,nStp);
z_t_hat = zeros(3,nStp);
z_err = zeros(3,nStp);

q = zeros(1,nP);
H_t = zeros(3,3,nStp);
S_t = zeros(3,3,nStp);
K_t = zeros(3,3,nStp);

Q_t = diag([app.sigma_d,app.sigma_phi,app.sigma_s]);

mu_bar_unc = zeros(3,1,nStp);
Sigma_bar_unc = zeros(3,3,nStp);

mu_bar_unc(:,:,1) = x_t_array(:,1);
Sigma_bar_unc(:,:,1) = zeros(3,3);

for k = 1:nStp-1
    u_t = v_t_array(:,k);
    mu_bar_unc(:,:,k+1) = g_t(app,u_t,mu_bar_unc(:,:,k));
    Sigma_bar_unc(:,:,k+1) = G_t(app,u_t,mu_bar_unc(:,:,k))*Sigma_bar_unc(:,:,k)*G_t(app,u_t,mu_bar_unc(:,:,k))'+V_t(app,u_t,mu_bar_unc(:,:,k))*M_t(app,u_t)*V_t(app,u_t,mu_bar_unc(:,:,k))';
end

% step 1
k = 1;

figure
grid on
hold on
axis equal
axis padded
for j = 1:nP
    plot(xM(j),yM(j),char(sM(j)))
end
plot(x_t_array(1,1),x_t_array(2,1),'.k')
plot(x_t_array(1,2),x_t_array(2,2),'.k')
plot(x_t_array(1,3),x_t_array(2,3),'.k')
drawEllipse(mu_bar_unc(:,:,1),Sigma_bar_unc(:,:,1),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,2),Sigma_bar_unc(:,:,2),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,3),Sigma_bar_unc(:,:,3),'#00FFFF','*','-')
drawRobot(x_t(:,k),0.1)
drawEllipse(mu_bar(:,:,k),Sigma_bar(:,:,k),'#000000','*','-')
drawEllipse(mu(:,:,k),Sigma(:,:,k),'#FF00FF','*','-')
axis([-1.25 1.25 -1.6 0.2])

% step 1 prediction
u_t = v_t_array(:,k);
mu_bar(:,:,k+1) = g_t(app,u_t,mu_bar(:,:,k));
Sigma_bar(:,:,k+1) = G_t(app,u_t,mu_bar(:,:,k))*Sigma_bar(:,:,k)*G_t(app,u_t,mu_bar(:,:,k))'+V_t(app,u_t,mu_bar(:,:,k))*M_t(app,u_t)*V_t(app,u_t,mu_bar(:,:,k))';

mu(:,:,k+1) = mu_bar(:,:,k+1);
Sigma(:,:,k+1) = Sigma_bar(:,:,k+1);

% step 2 measurement
k = 2;

m(:,k) = app.L_info(:,k-1);
z_t(:,k) = h_t(app,x_t(:,k),m(:,k));
m_hat(:,k) = [  x_t(1,k)+z_t(1,k)*cos(x_t(3,k)+z_t(2,k))
                x_t(2,k)+z_t(1,k)*sin(x_t(3,k)+z_t(2,k))
                z_t(3,k)];

z_t_hat(:,k) = h_t(app,mu_bar(:,:,k),m(:,k));
z_err(:,k) = z_t(:,k)-z_t_hat(:,k);

q(:,k) = (m(1,k)-mu_bar(1,:,k))^2+(m(2,k)-mu_bar(2,:,k))^2;
        
H_t(:,:,k) = [ -(m(1,k)-mu_bar(1,:,k))/sqrt(q(:,k)), -(m(2,k)-mu_bar(2,:,k))/sqrt(q(:,k)), 0
                (m(2,k)-mu_bar(2,:,k))/q(:,k), -(m(1,k)-mu_bar(1,:,k))/q(:,k), -1
                0, 0, 0];

S_t(:,:,k) = H_t(:,:,k)*Sigma_bar(:,:,k)*H_t(:,:,k)'+Q_t;
K_t(:,:,k) = Sigma_bar(:,:,k)*H_t(:,:,k)'/S_t(:,:,k);

mu(:,:,k) = mu_bar(:,:,k) + K_t(:,:,k)*z_err(:,k);
Sigma(:,:,k) = (eye(3)-K_t(:,:,k)*H_t(:,:,k))*Sigma_bar(:,:,k);

figure
grid on
hold on
axis equal
axis padded
for j = 1:nP
    plot(xM(j),yM(j),char(sM(j)))
end
plot(x_t_array(1,1),x_t_array(2,1),'.k')
plot(x_t_array(1,2),x_t_array(2,2),'.k')
plot(x_t_array(1,3),x_t_array(2,3),'.k')
drawEllipse(mu_bar_unc(:,:,1),Sigma_bar_unc(:,:,1),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,2),Sigma_bar_unc(:,:,2),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,3),Sigma_bar_unc(:,:,3),'#00FFFF','*','-')
drawRobot(x_t(:,k),0.1)
drawEllipse(mu_bar(:,:,k),Sigma_bar(:,:,k),'#000000','*',':')
drawEllipse(mu(:,:,k),Sigma(:,:,k),'#FF00FF','*','-')
plot([x_t(1,k) m_hat(1,k)],[x_t(2,k) m_hat(2,k)],'--r')
plot([mu_bar(1,:,k) m_hat(1,k)],[mu_bar(2,:,k) m_hat(2,k)],':r')
drawEllipse(m_hat(:,k),Q_t,'#FF0000','o','--')
axis([-1.25 1.25 -1.6 0.2])

% step 2 prediction
mu_bar(:,:,k) = mu(:,:,k);
Sigma_bar(:,:,k) = Sigma(:,:,k);

u_t = v_t_array(:,k);
mu_bar(:,:,k+1) = g_t(app,u_t,mu_bar(:,:,k));
Sigma_bar(:,:,k+1) = G_t(app,u_t,mu_bar(:,:,k))*Sigma_bar(:,:,k)*G_t(app,u_t,mu_bar(:,:,k))'+V_t(app,u_t,mu_bar(:,:,k))*M_t(app,u_t)*V_t(app,u_t,mu_bar(:,:,k))';

mu(:,:,k+1) = mu_bar(:,:,k+1);
Sigma(:,:,k+1) = Sigma_bar(:,:,k+1);


k = 3;

m(:,k) = app.L_info(:,k-1);
z_t(:,k) = h_t(app,x_t(:,k),m(:,k))
m_hat(:,k) = [  x_t(1,k)+z_t(1,k)*cos(x_t(3,k)+z_t(2,k))
                x_t(2,k)+z_t(1,k)*sin(x_t(3,k)+z_t(2,k))
                z_t(3,k)];

z_t_hat(:,k) = h_t(app,mu_bar(:,:,k),m(:,k))
z_err(:,k) = z_t(:,k)-z_t_hat(:,k);

q(:,k) = (m(1,k)-mu_bar(1,:,k))^2+(m(2,k)-mu_bar(2,:,k))^2;
        
H_t(:,:,k) = [ -(m(1,k)-mu_bar(1,:,k))/sqrt(q(:,k)), -(m(2,k)-mu_bar(2,:,k))/sqrt(q(:,k)), 0
                (m(2,k)-mu_bar(2,:,k))/q(:,k), -(m(1,k)-mu_bar(1,:,k))/q(:,k), -1
                0, 0, 0];

S_t(:,:,k) = H_t(:,:,k)*Sigma_bar(:,:,k)*H_t(:,:,k)'+Q_t;
K_t(:,:,k) = Sigma_bar(:,:,k)*H_t(:,:,k)'/S_t(:,:,k);

mu(:,:,k) = mu_bar(:,:,k) + K_t(:,:,k)*z_err(:,k);
Sigma(:,:,k) = (eye(3)-K_t(:,:,k)*H_t(:,:,k))*Sigma_bar(:,:,k);

figure
grid on
hold on
axis equal
axis padded
for j = 1:nP
    plot(xM(j),yM(j),char(sM(j)))
end
plot(x_t_array(1,1),x_t_array(2,1),'.k')
plot(x_t_array(1,2),x_t_array(2,2),'.k')
plot(x_t_array(1,3),x_t_array(2,3),'.k')
drawEllipse(mu_bar_unc(:,:,1),Sigma_bar_unc(:,:,1),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,2),Sigma_bar_unc(:,:,2),'#00FFFF','*','-')
drawEllipse(mu_bar_unc(:,:,3),Sigma_bar_unc(:,:,3),'#00FFFF','*','-')
drawRobot(x_t(:,k),0.1)
drawEllipse(mu_bar(:,:,k),Sigma_bar(:,:,k),'#000000','*',':')
drawEllipse(mu(:,:,k-1),Sigma(:,:,k-1),'#FF00FF','*','-')
drawEllipse(mu(:,:,k),Sigma(:,:,k),'#FF00FF','*','-')
plot([x_t(1,k) m_hat(1,k)],[x_t(2,k) m_hat(2,k)],'--r')
plot([mu_bar(1,:,k) m_hat(1,k)],[mu_bar(2,:,k) m_hat(2,k)],':r')
drawEllipse(m_hat(:,k),Q_t,'#FF0000','^','--')
axis([-1.25 1.25 -1.6 0.2])

% figure
% grid on
% hold on
% axis equal
% for k = 2:nStp-1
% %     drawEllipse(mu(:,k),Sigma(:,:,k),cmap(k, :))
% %     drawEllipse(mu_bar_e(:,k),Sigma_bar_e(:,:,k),'c')
%     drawEllipse(mu_bar(:,:,k),Sigma_bar(:,:,k),'k')
%     drawEllipse(mu_t_bar,Sigma_t_bar,'m')
% end

function drawEllipse(mu,Sigma,color,mrk,lin)

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

        plot(mu(1),mu(2),'Marker',mrk,'Color',color);
        plot(r_ellipse(:,1) + mu(1),r_ellipse(:,2) + mu(2),lin,'Color',color);
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

function out = h_t(app,x_t,m_j)

    tmp = [ sqrt((m_j(1)-x_t(1))^2+(m_j(2)-x_t(2))^2)
            atan2(m_j(2)-x_t(2),m_j(1)-x_t(1))-x_t(3)
            m_j(3) ];
    
    tmp = tmp + [sample(app,app.sigma_d^2);sample(app,app.sigma_phi^2);sample(app,app.sigma_s^2)];
    
    tmp(2) = normalizeFunc(tmp(2),-pi,pi);
    
    out = tmp;
end

function drawRobot(x_vec,r)

    xR = x_vec(1);
    yR = x_vec(2);
    tR = x_vec(3);

    plot(xR + r*cos(linspace(0,2*pi)),yR + r*sin(linspace(0,2*pi)),'-k',[xR,xR + r*cos(tR)],[yR,yR + r*sin(tR)],'-k',xR,yR,'.k')

end
        