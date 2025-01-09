clc, clear all, close all

app.nSmp = 1000;

app.dt = 1;
app.dist_type = 0;
app.prob_type = 0;
app.traj_type = 0;

app.alpha1 = 0.01;
app.alpha2 = 0.01;
app.alpha3 = 0.01;
app.alpha4 = 0.01;
app.alpha5 = 0.01;
app.alpha6 = 0.01;

% app.alpha1 = 0.0;
% app.alpha2 = 0.0;
% app.alpha3 = 0.0;
% app.alpha4 = 0.0;
% app.alpha5 = 0.0;
% app.alpha6 = 0.0;
    
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
        theta_array(i+1) = theta_array(i) + omega(i)*app.dt;
%         theta_array(i+1) = normalizeFunc(theta_array(i) + omega(i)*app.dt,-pi,pi);
%         if theta_array(i+1) == -pi
%             theta_array(i+1) = pi;
%         end
    end
else
    dt = 1;
    
    nStp = 2;

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

figure
plot(x_array,y_array,'.')
axis equal

x_grid = -2:0.01:2;
y_grid = -2:0.01:2;
[x_mesh,y_mesh] = meshgrid(x_grid,y_grid);

p_x_t = zeros(length(x_grid),length(y_grid),nStp-1);
for k = 1:nStp-1
    for i = 1:length(x_grid)
        for j = 1:length(y_grid)
            x_t = [x_mesh(i,j);y_mesh(i,j);x_t_array(3,k+1)];
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
drawRobot(x_t_array(:,1))
drawRobot(x_t_array(:,2))

for k = 1:nStp-1
    figure
    plot(x_t_array(1,:),x_t_array(2,:),'.k')
    hold on
    azSurf1 = surf(x_mesh,y_mesh,p_x_t(:,:,k));
    xlabel('x')
    ylabel('y')
    shading interp
    set(azSurf1, 'FaceAlpha', 'interp','AlphaData',p_x_t(:,:,k));
end

function drawRobot(x_vec)

    xR = x_vec(1);
    yR = x_vec(2);
    tR = x_vec(3);

    r = 0.1;

    hold("on")
    plot(xR + r*cos(linspace(0,2*pi)),yR + r*sin(linspace(0,2*pi)),'-k',[xR,xR + r*cos(tR)],[yR,yR + r*sin(tR)],'-k',xR,yR,'.k')

end

function n = normalizeFunc(value,start,finish) 

  width       = finish - start   ;
  offsetValue = value - start ;

  n = ( offsetValue - ( floor( offsetValue / width ) * width ) ) + start ;
  
end