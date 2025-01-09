
clc, clear all, close all

app.traj_type = 0;
app.prob_type = 0;
model = 1;

%% initialization
syms x y theta v omega real

syms delta_rot1 delta_trans delta_rot2

% state vector in symbolic form
x_vec = [x;y;theta];

% system funtions
g = [   x - (v/omega)*(sin(theta) - sin(theta+omega))
        y + (v/omega)*(cos(theta) - cos(theta+omega))
        theta + omega                                 ]   % nonlinear system model

% delta_rot1  = atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar;
% delta_trans = sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2);
% delta_rot2  = theta_bar_prm - theta_bar - delta_rot1;
    
    x_prm       = x + delta_trans*cos(theta + delta_rot1);
    y_prm       = y + delta_trans*sin(theta + delta_rot1);
    theta_prm   = theta + delta_rot1 + delta_rot2;
    
% system funtions
g2 = [	x + delta_trans*cos(theta + delta_rot1)
        y + delta_trans*sin(theta + delta_rot1)
        theta + delta_rot1 + delta_rot2                                ]   % nonlinear system model

% linearization of the functions
G = [diff(g,x) diff(g,y) diff(g,theta)]

% linearization of the functions
G2 = [diff(g2,x) diff(g2,y) diff(g2,theta)]

gb = [ 	x + v*cos(theta + omega)
        y + v*sin(theta + omega)
        theta + omega                	]   % nonlinear system model

% linearization of the functions
Gb = [diff(gb,x) diff(gb,y) diff(gb,theta)]

dt = 1;

%% generate trajectory
if app.traj_type == 0
    nStp = ceil(8/dt);

    x = 1;
    y = 0;
    theta = -pi/2;

    v = pi/4*ones(1,nStp); % m/s
    omega = -pi/4*ones(1,nStp); % rad/s

    x_array     = x*ones(1,nStp);
    y_array     = y*ones(1,nStp);
    theta_array = theta*ones(1,nStp);

    for i = 1:nStp-1
        x_array(i+1) = x_array(i) - (v(i)/omega(i))*(sin(theta_array(i)) - sin(theta_array(i)+omega(i)*dt));
        y_array(i+1) = y_array(i) + (v(i)/omega(i))*(cos(theta_array(i)) - cos(theta_array(i)+omega(i)*dt));
%         x_array(i+1) = x_array(i) + v(i)*dt*(cos(theta_array(i)+omega(i)*dt));
%         y_array(i+1) = y_array(i) + v(i)*dt*(sin(theta_array(i)+omega(i)*dt));
        theta_array(i+1) = theta_array(i) + omega(i)*dt;
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

figure
plot(x_array,y_array,'.')
axis equal


% system funtions
g_t = @(x,y,theta,v,omega)[ x - (v/omega)*(sin(theta) - sin(theta+omega))
                            y + (v/omega)*(cos(theta) - cos(theta+omega))
                            theta + omega                                   ];   % nonlinear system model
                        
G_t = @(theta,v,omega)[ 1,0,(v/omega)*(cos(omega + theta) - cos(theta))
                        0,1,(v/omega)*(sin(omega + theta) - sin(theta))
                        0,0,1                                               ];   % model jacobian

g_t1 = g_t(x_array(1),y_array(1),theta_array(1),v_t_array(1,1),v_t_array(2,1))
G_t1 = G_t(theta_array(1),v_t_array(1,1),v_t_array(2,1))

g2_t = @(x,y,theta,x_bar,y_bar,theta_bar,x_bar_prm,y_bar_prm,theta_bar_prm)[	x + sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*cos(theta + atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar)
                                                                                y + sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*sin(theta + atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar)
                                                                                theta + theta_bar_prm - theta_bar       ];   % nonlinear system model

G2_t = @(theta,x_bar,y_bar,theta_bar,x_bar_prm,y_bar_prm)[	1, 0, -sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*sin(atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar + theta)
                                                            0, 1,  sqrt((x_bar-x_bar_prm)^2 + (y_bar-y_bar_prm)^2)*cos(atan2(y_bar_prm-y_bar, x_bar_prm-x_bar) - theta_bar + theta)
                                                            0, 0,                                    1];
 
g2_t1 = g2_t(x_array(1),y_array(1),theta_array(1),x_array(1),y_array(1),theta_array(1),x_array(2),y_array(2),theta_array(2))
G2_t1 = G2_t(theta_array(1),x_array(1),y_array(1),theta_array(1),x_array(2),y_array(2))
