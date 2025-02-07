clc, close all, clear all


nStp = 13;

x_array = [0,2,4,6,8,8,8,8,8,6,4,2,0]/2;
y_array = [0,0,0,0,0,2,4,6,8,8,8,8,8]/2;
theta_array = [0,0,0,0,0,pi/2,pi/2,pi/2,pi,pi,pi,pi,pi];

x_r = x_array(1,1);
y_r = y_array(1,1);
theta_r = theta_array(1,1);

v = zeros(1,nStp); % m/s
omega = zeros(1,nStp); % rad/s

for i = 1:nStp-1
    v(1,i) = sqrt((x_array(1,i+1)-x_array(1,i))^2+(y_array(1,i+1)-y_array(1,i))^2); % m/s
    omega(1,i) = theta_array(1,i+1)-theta_array(1,i); % rad/s
end

x_t_array = [x_array;y_array;theta_array];
v_t_array = [v;omega];

dt = 0.5;
a = 0;
k = 1;
x_s = zeros(1,floor(nStp/dt)-1);
y_s = zeros(1,floor(nStp/dt)-1);
theta_s = zeros(1,floor(nStp/dt)-1);
for i = 1:nStp
    
    if i < nStp
        dx = (x_array(i+1)-x_array(i))*dt;
        dy = (y_array(i+1)-y_array(i))*dt;
        dtheta = (theta_array(i+1)-theta_array(i))*dt;
    else
        a = 1;
        dx = (x_array(i)-x_array(i))*dt;
        dy = (y_array(i)-y_array(i))*dt;
        dtheta = (theta_array(i)-theta_array(i))*dt;
    end
    
    for j = 1:floor(1/dt)-a
        
        x_s(k) = x_array(i)+(j-1)*dx;
        y_s(k) = y_array(i)+(j-1)*dy;
        theta_s(k) = theta_array(i)+(j-1)*dtheta;
        
        k = k+1;
    end
    
end
% x_s(k+1) = x_array(i+1);
% y_s(k+1) = y_array(i+1);
% theta_s(k+1) = theta_array(i+1);

% for j = 1:floor(1/dt)
% 
%     x_s(k) = x_array(i)+(j-1)*dx;
%     y_s(k) = y_array(i)+(j-1)*dy;
%     theta_s(k) = theta_array(i)+(j-1)*dtheta;
% 
%     k = k+1;
% end
x_s

figure
hold on
for k = 1:floor(nStp/dt)-1
    plot(x_s(k),y_s(k),'.')
    pause(0.5)
end