% KOM613E Homework #2a
% Mertcan Kaya
% 518192004

clc, clear all, close all

%% initialization
format long
syms x y theta real

% state vector in symbolic form
x_vec = [x y theta];

% initial estimate
mu_t0 = [0 0 0];
Sigma_t0 = [0.01 0 0;0 1000 0;0 0 0.1];

% disturbances
R = zeros(3); % no model noise
Q = 0.01;   % measurement noise

% system funtions
g = [x+cos(theta);y+sin(theta);theta]   % nonlinear system model
h = x % measurement model

% linearization of the functions
G = [diff(g,x) diff(g,y) diff(g,theta)]
H = [diff(h,x) diff(h,y) diff(h,theta)]

%% first time step

% possible initial states
x_vec_t0 = mvnrnd(mu_t0,Sigma_t0,100)';
x_t0 = x_vec_t0(1,:);
y_t0 = x_vec_t0(2,:);
theta_t0 = x_vec_t0(3,:);

%% second time step

% possible states
posNb = 100;
ornNb = 10;
x_vec_t1_tmp = zeros(3,posNb,ornNb);
for i = 1:posNb
    for j = 1:ornNb
        x_vec_t1_tmp(:,i,j) = double(subs(g,x_vec,[x_t0(i) y_t0(i) theta_t0(j)]));
    end
end

x_vec_t1_pst = reshape(x_vec_t1_tmp,3,posNb*ornNb);

H_t0 = double(subs(H,x_vec,mu_t0))

z_t1_pst = H_t0*x_vec_t1_pst;

% select a random possible state as the actual state
rand_pnt = floor(posNb*ornNb*rand)+1;
delta = normrnd(0,Q);
x_t1 = x_vec_t1_pst(:,rand_pnt)
z_t1 = z_t1_pst(rand_pnt)+delta

%% Extended Kalman Filter

g_t0 = double(subs(g,x_vec,mu_t0))
G_t0 = double(subs(G,x_vec,mu_t0))

% prediction
mu_t1_bar       = g_t0
Sigma_t1_bar    = G_t0*Sigma_t0*G_t0'+R

I = eye(size(Sigma_t1_bar));

% correction
K_t1            = Sigma_t1_bar*H_t0'/(H_t0*Sigma_t1_bar*H_t0'+Q)
Sigma_t1        = (I-K_t1*H_t0)*Sigma_t1_bar
mu_t1           = mu_t1_bar + K_t1*(z_t1-H_t0*mu_t1_bar)

%% calculate first time step

% the 95% confidence interval error ellipse
chisquare_val = sqrt(5.991);
theta_grid = linspace(0,2*pi);

[eigenvec_t0, eigenval_t0] = eig(Sigma_t0(1:2,1:2));

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval_t0 == max(max(eigenval_t0)));
largest_eigenvec_t0 = eigenvec_t0(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t0 = max(max(eigenval_t0));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t0 = max(eigenval_t0(:,2));
    smallest_eigenvec_t0 = eigenvec_t0(:,2);
else
    smallest_eigenval_t0 = max(eigenval_t0(:,1));
    smallest_eigenvec_t0 = eigenvec_t0(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t0 = atan2(largest_eigenvec_t0(2), largest_eigenvec_t0(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t0 < 0)
    angle_t0 = angle_t0 + 2*pi;
end

phi_t0 = angle_t0;
a_t0 = chisquare_val * sqrt(largest_eigenval_t0);
b_t0 = chisquare_val * sqrt(smallest_eigenval_t0);

% the ellipse in x and y coordinates 
ellipse_x_r_t0  = a_t0*cos(theta_grid);
ellipse_y_r_t0  = b_t0*sin(theta_grid);

% rotation matrix of angle phi
rot_t0 = [ cos(phi_t0) sin(phi_t0); -sin(phi_t0) cos(phi_t0) ];

% rotation of the ellipse to angle phi
r_ellipse_t0 = [ellipse_x_r_t0;ellipse_y_r_t0]' * rot_t0;

%% calculate possible

% uncertainty ellipse of possible states
mu_t1_psb = mean(x_vec_t1_pst(1:2,:),2);
sigma_t1_psb = std(x_vec_t1_pst(1:2,:),1,2); % standard deviation
Sigma_t1_psb = diag(sigma_t1_psb.^2); % covariance

[eigenvec_t1_psb, eigenval_t1_psb] = eig(Sigma_t1_psb(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval_t1_psb == max(max(eigenval_t1_psb)));
largest_eigenvec_t1_psb = eigenvec_t1_psb(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1_psb = max(max(eigenval_t1_psb));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1_psb = max(eigenval_t1_psb(:,2));
    smallest_eigenvec_t1_psb = eigenvec_t1_psb(:,2);
else
    smallest_eigenval_t1_psb = max(eigenval_t1_psb(:,1));
    smallest_eigenvec_t1_psb = eigenvec_t1_psb(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1_psb = atan2(largest_eigenvec_t1_psb(2), largest_eigenvec_t1_psb(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1_psb < 0)
    angle_t1_psb = angle_t1_psb + 2*pi;
end

phi_t1_psb = angle_t1_psb;
a_t1_psb = chisquare_val * sqrt(largest_eigenval_t1_psb);
b_t1_psb = chisquare_val * sqrt(smallest_eigenval_t1_psb);

% the ellipse in x and y coordinates 
ellipse_x_r_t1_psb  = a_t1_psb*cos(theta_grid);
ellipse_y_r_t1_psb  = b_t1_psb*sin(theta_grid);

% rotation matrix of angle phi
rot_t1_psb = [cos(phi_t1_psb) sin(phi_t1_psb); -sin(phi_t1_psb) cos(phi_t1_psb)];

% rotation of the ellipse to angle phi
r_ellipse_t1_psb = [ellipse_x_r_t1_psb;ellipse_y_r_t1_psb]' * rot_t1_psb;

%% calculate predicted

% draw predicted states in the second time step
x_vec_t1_bar = mvnrnd(mu_t1_bar,Sigma_t1_bar,100)';

% calculate the uncertainty ellipse of predicted states
[eigenvec_t1_bar, eigenval_t1_bar] = eig(Sigma_t1_bar(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval_t1_bar == max(max(eigenval_t1_bar)));
largest_eigenvec_t1_bar = eigenvec_t1_bar(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1_bar = max(max(eigenval_t1_bar));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1_bar = max(eigenval_t1_bar(:,2));
    smallest_eigenvec_t1_bar = eigenvec_t1_bar(:,2);
else
    smallest_eigenval_t1_bar = max(eigenval_t1_bar(:,1));
    smallest_eigenvec_t1_bar = eigenvec_t1_bar(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1_bar = atan2(largest_eigenvec_t1_bar(2), largest_eigenvec_t1_bar(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1_bar < 0)
    angle_t1_bar = angle_t1_bar + 2*pi;
end

phi_t1_bar = angle_t1_bar;
a_t1_bar = chisquare_val * sqrt(largest_eigenval_t1_bar);
b_t1_bar = chisquare_val * sqrt(smallest_eigenval_t1_bar);

% the ellipse in x and y coordinates 
ellipse_x_r_t1_bar  = a_t1_bar * cos(theta_grid);
ellipse_y_r_t1_bar  = b_t1_bar * sin(theta_grid);

% rotation matrix of angle phi
rot_t1_bar = [cos(phi_t1_bar) sin(phi_t1_bar); -sin(phi_t1_bar) cos(phi_t1_bar)];

% rotation of the ellipse to angle phi
r_ellipse_t1_bar = [ellipse_x_r_t1_bar;ellipse_y_r_t1_bar]' * rot_t1_bar;

x_vec_t1_bar_zoom = mvnrnd(mu_t1_bar,Sigma_t1_bar,100)';

%% calculate corrected

[eigenvec_t1, eigenval_t1] = eig(Sigma_t1(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval_t1 == max(max(eigenval_t1)));
largest_eigenvec_t1 = eigenvec_t1(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1 = max(max(eigenval_t1));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1 = max(eigenval_t1(:,2));
    smallest_eigenvec_t1 = eigenvec_t1(:,2);
else
    smallest_eigenval_t1 = max(eigenval_t1(:,1));
    smallest_eigenvec_t1 = eigenvec_t1(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1 = atan2(largest_eigenvec_t1(2), largest_eigenvec_t1(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1 < 0)
    angle_t1 = angle_t1 + 2*pi;
end

phi_t1 = angle_t1;
a_t1 = chisquare_val * sqrt(largest_eigenval_t1);
b_t1 = chisquare_val * sqrt(smallest_eigenval_t1);

% the ellipse in x and y coordinates 
ellipse_x_r_t1 = a_t1 * cos(theta_grid);
ellipse_y_r_t1 = b_t1 * sin(theta_grid);

% rotation matrix of angle phi
rot_t1 = [cos(phi_t1) sin(phi_t1); -sin(phi_t1) cos(phi_t1)];

% rotation of the ellipse to angle phi
r_ellipse_t1 = [ellipse_x_r_t1;ellipse_y_r_t1]' * rot_t1;

%% draw gaussian probability distribution projected on the x-axis

x_grid = linspace(-0.5,1.5,400);
x_vec_grid = zeros(3,length(x_grid));
x_vec_grid(1,:) = x_grid;

% p_x_t0_vec = zeros(3,length(x_grid));
p_x_t0 = zeros(1,length(x_grid));
% p_x_t1_vec_bar = zeros(3,length(x_grid));
p_x_t1_psb = zeros(1,length(x_grid));
p_x_t1_bar = zeros(1,length(x_grid));
p_z_t1 = zeros(1,length(x_grid));
% p_x_t1_vec = zeros(3,length(x_grid));
p_x_t1 = zeros(1,length(x_grid));
for i = 1:length(x_grid)
%     p_x_t0_vec(:,i) = det(2*pi*sigma_t0)^(-1/2)*exp((-1/2)*(x_vec_grid(:,i)-mu_t0)'*(sigma_t0\(x_vec_grid(:,i)-mu_t0)));
    p_x_t0(:,i) = (2*pi*Sigma_t0(1,1))^(-1/2)*exp((-1/2)*(x_grid(:,i)-mu_t0(1))^2/Sigma_t0(1,1));
%     p_x_t1_vec_bar(:,i) = det(2*pi*sigma_t1_bar)^(-1/2)*exp((-1/2)*(x_vec_grid(:,i)-mu_t1_bar)'*(sigma_t1_bar\(x_vec_grid(:,i)-mu_t1_bar)));
    p_x_t1_psb(:,i) = (2*pi*Sigma_t1_psb(1,1))^(-1/2)*exp((-1/2)*(x_grid(:,i)-mu_t1_psb(1))^2/Sigma_t1_psb(1,1));
    p_x_t1_bar(:,i) = (2*pi*Sigma_t1_bar(1,1))^(-1/2)*exp((-1/2)*(x_grid(:,i)-mu_t1_bar(1))^2/Sigma_t1_bar(1,1));
    p_z_t1(:,i) = (2*pi*Q)^(-1/2)*exp((-1/2)*(x_grid(:,i)-z_t1)^2/Q);
%     p_x_t1_vec(:,i) = det(2*pi*sigma_t1)^(-1/2)*exp((-1/2)*(x_vec_grid(:,i)-mu_t1)'*(sigma_t1\(x_vec_grid(:,i)-mu_t1)));
    p_x_t1(:,i) = (2*pi*Sigma_t1(1,1))^(-1/2)*exp((-1/2)*(x_grid(:,i)-mu_t1(1))^2/Sigma_t1(1,1));
end

%% figure 1

fig1 = figure(1);
ax1 = axes('Parent',fig1);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-0.5, 1.5]);
ylim([-100, 100]);

plt_t0(1) = plot(x_t0,y_t0,'.b');
plt_t0(2) = plot(r_ellipse_t0(:,1) + mu_t0(1),r_ellipse_t0(:,2) + mu_t0(2),'-b');
plt_t0(3) = plot(mu_t0(1),mu_t0(2),'*b','LineWidth',2);

plt_t1_pst(1) = plot(x_vec_t1_pst(1,:),x_vec_t1_pst(2,:),'.g');
plt_t1_pst(2) = plot(r_ellipse_t1_psb(:,1) + mu_t1_psb(1),r_ellipse_t1_psb(:,2) + mu_t1_psb(2),'-g');
plt_t1_pst(3) = plot(mu_t1_psb(1),mu_t1_psb(2),'*g','LineWidth',2);

plt_t1_bar(1) = plot(x_vec_t1_bar_zoom(1,:),x_vec_t1_bar_zoom(2,:),'.r');
plt_t1_bar(2) = plot(r_ellipse_t1_bar(:,1) + mu_t1_bar(1),r_ellipse_t1_bar(:,2) + mu_t1_bar(2),'-r');
plt_t1_bar(3) = plot(mu_t1_bar(1),mu_t1_bar(2),'*r','LineWidth',2);

plt_t1_act = plot(x_t1(1),x_t1(2),'*k','LineWidth',2);

% h = legend([plt_t0,plt_t1_pst,plt_t1_bar,plt_t1_act],'$\overline{bel}(x_1)$','$\overline{UE}_{1,95\%}$','$\bar{\mu}_1$','$x_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 2

% fig2 = figure(2);
% ax2 = axes('Parent',fig2);
% hold on
% grid on
% xlabel('x')
% ylabel('y')
% 
% copyobj(p_t0,ax2)
% copyobj(e_t0,ax2)
% copyobj(p_t1_psb,ax2)
% copyobj(e_t1_psb,ax2)
% copyobj(p_t1_slc,ax2)
% plot(ax2,x_vec_t1_bar_zoom(1,:),x_vec_t1_bar_zoom(2,:),'r.');
% copyobj(e_t1_bar,ax2)

%% figure 3

fig3 = figure(3);
ax3 = axes('Parent',fig3);
hold on
grid on
xlabel('x')
ylabel('y')
xlim([-0.5, 1.5]);
ylim([-100, 100]);

% plot the eigenvectors
plt_qvr(1) = quiver(mu_t0(1), mu_t0(2), eigenvec_t0(1,1)*sqrt(eigenval_t0(1,1)), eigenvec_t0(2,1)*sqrt(eigenval_t0(1,1)), '-b', 'LineWidth', 2, 'ShowArrowHead','off');
plt_qvr(2) = quiver(mu_t0(1), mu_t0(2), eigenvec_t0(1,2)*sqrt(eigenval_t0(2,2)), eigenvec_t0(2,2)*sqrt(eigenval_t0(2,2)), '-b', 'LineWidth', 2, 'ShowArrowHead','off');

% plot the eigenvectors
plt_qvr(3) = quiver(mu_t1_psb(1), mu_t1_psb(2), eigenvec_t1_psb(1,1)*sqrt(eigenval_t1_psb(1,1)), eigenvec_t1_psb(2,1)*sqrt(eigenval_t1_psb(1,1)), '-g', 'LineWidth', 2, 'ShowArrowHead','off');
plt_qvr(4) = quiver(mu_t1_psb(1), mu_t1_psb(2), eigenvec_t1_psb(1,2)*sqrt(eigenval_t1_psb(2,2)), eigenvec_t1_psb(2,2)*sqrt(eigenval_t1_psb(2,2)), '-g', 'LineWidth', 2, 'ShowArrowHead','off');

% plot the eigenvectors
plt_qvr(5) = quiver(mu_t1_bar(1), mu_t1_bar(2), largest_eigenvec_t1_bar(1)*sqrt(largest_eigenval_t1_bar), largest_eigenvec_t1_bar(2)*sqrt(largest_eigenval_t1_bar), '-r', 'LineWidth', 2, 'ShowArrowHead', 'off');
plt_qvr(6) = quiver(mu_t1_bar(1), mu_t1_bar(2), smallest_eigenvec_t1_bar(1)*sqrt(smallest_eigenval_t1_bar), smallest_eigenvec_t1_bar(2)*sqrt(smallest_eigenval_t1_bar), '-r', 'LineWidth', 2, 'ShowArrowHead', 'off');

% plot the eigenvectors
plt_qvr(7) = quiver(mu_t1(1), mu_t1(2), largest_eigenvec_t1(1)*sqrt(largest_eigenval_t1), largest_eigenvec_t1(2)*sqrt(largest_eigenval_t1), '-m', 'LineWidth', 2, 'ShowArrowHead', 'off');
plt_qvr(8) = quiver(mu_t1(1), mu_t1(2), smallest_eigenvec_t1(1)*sqrt(smallest_eigenval_t1), smallest_eigenvec_t1(2)*sqrt(smallest_eigenval_t1), '-m', 'LineWidth', 2, 'ShowArrowHead', 'off');

copyobj(plt_t0(2:3),ax3)
copyobj(plt_t1_pst(2:3),ax3)
copyobj(plt_t1_bar(2:3),ax3)
copyobj(plt_t1_act,ax3)

plt_z_t1(1) = xline(z_t1,'--c');
plt_z_t1(2) = xline(z_t1-chisquare_val*sqrt(Q),'c');
plt_z_t1(3) = xline(z_t1+chisquare_val*sqrt(Q),'c');

% plot the error ellipse
plt_t1(1) = plot(r_ellipse_t1(:,1) + mu_t1(1),r_ellipse_t1(:,2) + mu_t1(2),'-m');
plt_t1(2) = plot(mu_t1(1),mu_t1(2),'*m','LineWidth',2);

%% figure 4

% fig4 = figure(4);
% ax4 = axes('Parent',fig4);
% hold on
% grid on
% xlabel('x')
% ylabel('y')
% 
% copyobj(e_t0,ax4)
% copyobj(e_t1_psb,ax4)
% copyobj(e_t1_bar,ax4)
% xline(z_t1-chisquare_val*sqrt(Q))
% xline(z_t1+chisquare_val*sqrt(Q))
% copyobj(e_t1,ax4)
% 
% plot(z_t1,0,'k*','LineWidth',2)
% plot(mu_t1_psb(1),0,'g*','LineWidth',2)
% plot(mu_t1_bar(1),0,'r*','LineWidth',2)
% plot(mu_t0(1),0,'b*','LineWidth',2)
% plot(mu_t1(1),0,'m*','LineWidth',2)

%% figure 5

fig4 = figure(5);
hold on
grid on
xlabel('x')
ylabel('p(x)')

plt_p(1) = plot(x_grid,p_x_t0,'b');
plot([mu_t0(1),mu_t0(1)],[0,(2*pi*Sigma_t0(1,1))^(-1/2)],'--b')
plt_p(2) = plot(x_grid,p_x_t1_psb,'g');
plot([mu_t1_psb(1),mu_t1_psb(1)],[0,(2*pi*Sigma_t1_psb(1,1))^(-1/2)],'--g')
plt_p(3) = plot(x_grid,p_x_t1_bar,'r');
plot([mu_t1_bar(1),mu_t1_bar(1)],[0,(2*pi*Sigma_t1_bar(1,1))^(-1/2)],'--r')
xline(x_t1(1),'--k');
plt_p(4) = plot(x_grid,p_z_t1,'c');
plot([z_t1,z_t1],[0,(2*pi*Q)^(-1/2)],'--c')
plt_p(5) = plot(x_grid,p_x_t1,'m');
plot([mu_t1(1),mu_t1(1)],[0,(2*pi*Sigma_t1(1,1))^(-1/2)],'--m')
