% KOM613E Homework #2a
% Mertcan Kaya
% 518192004

clc, clear all, close all

%% initialization
format long
syms x y theta real

% state vector in symbolic form
x_vec = [x;y;theta];

% initial estimate
mu_t0 = [0;0;0];
Sigma_t0 = [0.01 0 0;0 0.01 0;0 0 10000];

% disturbances
R = zeros(3); % no model noise
Q = 0.01;   % measurement noise

% system funtions
g = [x+cos(theta);y+sin(theta);theta]   % nonlinear system model
h = x % measurement model

% linearization of the functions
G = [diff(g,x) diff(g,y) diff(g,theta)]
H = [diff(h,x) diff(h,y) diff(h,theta)]

%% Sample points

% possible initial time step states
x_vec_t0 = mvnrnd(mu_t0,Sigma_t0,100)';
x_t0 = x_vec_t0(1,:);
y_t0 = x_vec_t0(2,:);
theta_t0 = x_vec_t0(3,:);

% possible second time step states
posNb = 20;
ornNb = 50;
x_vec_t1_tmp = zeros(3,posNb,ornNb);
for i = 1:posNb
    for j = 1:ornNb
        x_vec_t1_tmp(:,i,j) = double(subs(g,x_vec,[x_t0(i);y_t0(i);theta_t0(j)]));
    end
end
x_vec_t1 = reshape(x_vec_t1_tmp,3,posNb*ornNb);

% select a random possible state as the actual state
rand_pnt = floor(posNb*ornNb*rand)+1;
delta = normrnd(0,sqrt(Q));
x_t1 = x_vec_t1(:,rand_pnt)

% measurement
H_t0 = double(subs(H,x_vec,mu_t0))
z_t1 = H_t0*x_vec_t1(:,rand_pnt)+delta

% z_t1_vec = zeros(1,100);
% for i = 1:100
%     z_t1_vec(i) = H_t0*x_vec_t1(:,rand_pnt)+normrnd(0,sqrt(Q));
% end

%% Kalman Filter for Intuitive Solution

% prediction for intuitive posterior
mu_t1_mc    = mean(x_vec_t1,2);
sigma_t1_mc = std(x_vec_t1,1,2); % standard deviation
Sigma_t1_mc = diag(sigma_t1_mc.^2); % covariance

I = eye(size(Sigma_t1_mc));

% correction for intuitive posterior
K_t1_kf    	= Sigma_t1_mc*H_t0'/(H_t0*Sigma_t1_mc*H_t0'+Q)
Sigma_t1_kf	= (I-K_t1_kf*H_t0)*Sigma_t1_mc
mu_t1_kf   	= mu_t1_mc + K_t1_kf*(z_t1-H_t0*mu_t1_mc)

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

%% Unscented Kalman Filter

n = length(x_vec);
chiNb = 2*n + 1;

kappa = 10;      % kappa >= 0
alpha = 0.5;    % alpha € (0,1]
lambda = alpha^2*(n+kappa)-n;
beta = 2;

chi_t0 = zeros(3,chiNb);
mu_t0_offset = sqrtm((n+lambda)*Sigma_t0);
% gamma = sqrt(n+lambda)
% mu_offset2 = gamma*sqrtm(Sigma_t0)
% [V,D] = eig((n+lambda)*Sigma_t0)
% S = V*diag(sqrt(diag(D)))*V'
% chol((n+lambda)*Sigma_t0)
% gamma*chol(Sigma_t0)
for k = 1:chiNb
    if k == 1
        chi_t0(:,k) = mu_t0;
    elseif k <= n+1
        chi_t0(:,k) = mu_t0 + mu_t0_offset(:,k-1);
    else
        chi_t0(:,k) = mu_t0 - mu_t0_offset(:,k-n-1);
    end
end
chi_t0;

omega_m = zeros(1,chiNb);
omega_c = zeros(1,chiNb);
for k = 1:chiNb
    if k == 1
        omega_m(k) = lambda/(n+lambda);
        omega_c(k) = lambda/(n+lambda)+1-alpha^2+beta;
    else
        omega_m(k) = 1/(2*(n+lambda));
        omega_c(k) = 1/(2*(n+lambda));
    end
end

chi_t1_bar_s = zeros(3,chiNb);
for k = 1:chiNb
    chi_t1_bar_s(:,k) = double(subs(g,x_vec,chi_t0(:,k)));
end

mu_t1_bar_u = zeros(3,1);
for k = 1:chiNb
    mu_t1_bar_u = mu_t1_bar_u + omega_m(k)*chi_t1_bar_s(:,k);
end

Sigma_t1_bar_u = zeros(3,3);
for k = 1:chiNb
    Sigma_t1_bar_u = Sigma_t1_bar_u + omega_c(k)*(chi_t1_bar_s(:,k)-mu_t1_bar_u)*(chi_t1_bar_s(:,k)-mu_t1_bar_u)';
end
Sigma_t1_bar_u = Sigma_t1_bar_u + R;

chi_t1_bar = zeros(3,chiNb);
mu_t1_offset_bar = sqrtm((n+lambda)*Sigma_t1_bar_u);
for k = 1:chiNb
    if k == 1
        chi_t1_bar(:,k) = mu_t1_bar;
    elseif k <= n+1
        chi_t1_bar(:,k) = mu_t1_bar + mu_t1_offset_bar(:,k-1);
    else
        chi_t1_bar(:,k) = mu_t1_bar - mu_t1_offset_bar(:,k-n-1);
    end
end
chi_t1_bar;

zeta_t1_bar = zeros(1,chiNb);
for k = 1:chiNb
    zeta_t1_bar(:,k) = double(subs(h,x_vec,chi_t1_bar(:,k)));
end

z_t1_hat = 0;
for k = 1:chiNb
    z_t1_hat = z_t1_hat + omega_m(k)*zeta_t1_bar(:,k);
end

S_t1 = 0;
for k = 1:chiNb
    S_t1 = S_t1 + omega_c(k)*(zeta_t1_bar(:,k)-z_t1_hat)*(zeta_t1_bar(:,k)-z_t1_hat)';
end
S_t1 = S_t1 + Q;

Sigma_t1_bar_xz = zeros(3,1);
for k = 1:chiNb
    Sigma_t1_bar_xz = Sigma_t1_bar_xz + omega_c(k)*((chi_t1_bar(:,k)-mu_t1_bar_u)*(zeta_t1_bar(:,k)-z_t1_hat)');
end

K_t1_u = Sigma_t1_bar_xz/S_t1;

mu_t1_u = mu_t1_bar_u + K_t1_u*(z_t1-z_t1_hat)';
Sigma_t1_u = Sigma_t1_bar_u - K_t1_u*S_t1*K_t1_u';

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

%% calculate Monte-Carlo

[eigenvec_t1_psb, eigenval_t1_psb] = eig(Sigma_t1_mc(1:2,1:2))

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

%% calculate corrected monte-carlo

[eigenvec_t1_kf, eigenval_t1_kf] = eig(Sigma_t1_kf(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval_t1_kf == max(max(eigenval_t1_kf)));
largest_eigenvec_t1_kf = eigenvec_t1_kf(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1_kf = max(max(eigenval_t1_kf));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1_kf = max(eigenval_t1_kf(:,2));
    smallest_eigenvec_t1_kf = eigenvec_t1_kf(:,2);
else
    smallest_eigenval_t1_kf = max(eigenval_t1_kf(:,1));
    smallest_eigenvec_t1_kf = eigenvec_t1_kf(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1_kf = atan2(largest_eigenvec_t1_kf(2), largest_eigenvec_t1_kf(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1_kf < 0)
    angle_t1_kf = angle_t1_kf + 2*pi;
end

phi_t1_kf = angle_t1_kf;
a_t1_kf = chisquare_val * sqrt(largest_eigenval_t1_kf);
b_t1_kf = chisquare_val * sqrt(smallest_eigenval_t1_kf);

% the ellipse in x and y coordinates 
ellipse_x_r_t1_kf = a_t1_kf * cos(theta_grid);
ellipse_y_r_t1_kf = b_t1_kf * sin(theta_grid);

% rotation matrix of angle phi
rot_t1_kf = [cos(phi_t1_kf) sin(phi_t1_kf); -sin(phi_t1_kf) cos(phi_t1_kf)];

% rotation of the ellipse to angle phi
r_ellipse_t1_kf = [ellipse_x_r_t1_kf;ellipse_y_r_t1_kf]' * rot_t1_kf;

%% calculate predicted

% draw predicted states in the second time step
x_vec_t1_bar = mvnrnd(mu_t1_bar,Sigma_t1_bar,200)';

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

x_vec_t1_bar_zoom = mvnrnd(mu_t1_bar,Sigma_t1_bar,10000)';

%% calculate corrected for EKF

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

%% calculate predicted unscented transformation UE

[eigenvec_t1_bar_u, eigenval_t1_bar_u] = eig(Sigma_t1_bar_u(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval_t1_bar_u == max(max(eigenval_t1_bar_u)));
largest_eigenvec_t1_bar_u = eigenvec_t1_bar_u(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1_bar_u = max(max(eigenval_t1_bar_u));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1_bar_u = max(eigenval_t1_bar_u(:,2));
    smallest_eigenvec_t1_bar_u = eigenvec_t1_bar_u(:,2);
else
    smallest_eigenval_t1_bar_u = max(eigenval_t1_bar_u(:,1));
    smallest_eigenvec_t1_bar_u = eigenvec_t1_bar_u(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1_bar_u = atan2(largest_eigenvec_t1_bar_u(2), largest_eigenvec_t1_bar_u(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1_bar_u < 0)
    angle_t1_bar_u = angle_t1_bar_u + 2*pi;
end

phi_t1_bar_u = angle_t1_bar_u;
a_t1_bar_u = chisquare_val * sqrt(largest_eigenval_t1_bar_u);
b_t1_bar_u = chisquare_val * sqrt(smallest_eigenval_t1_bar_u);

% the ellipse in x and y coordinates 
ellipse_x_r_t1_bar_u  = a_t1_bar_u*cos(theta_grid);
ellipse_y_r_t1_bar_u  = b_t1_bar_u*sin(theta_grid);

% rotation matrix of angle phi
rot_t1_bar_u = [cos(phi_t1_bar_u) sin(phi_t1_bar_u); -sin(phi_t1_bar_u) cos(phi_t1_bar_u)];

% rotation of the ellipse to angle phi
r_ellipse_t1_bar_u = [ellipse_x_r_t1_bar_u;ellipse_y_r_t1_bar_u]' * rot_t1_bar_u;

%% calculate corrected unscented transformation UE

[eigenvec_t1_u, eigenval_t1_u] = eig(Sigma_t1_u(1:2,1:2))

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval_t1_u == max(max(eigenval_t1_u)));
largest_eigenvec_t1_u = eigenvec_t1_u(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval_t1_u = max(max(eigenval_t1_u));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval_t1_u = max(eigenval_t1_u(:,2));
    smallest_eigenvec_t1_u = eigenvec_t1_u(:,2);
else
    smallest_eigenval_t1_u = max(eigenval_t1_u(:,1));
    smallest_eigenvec_t1_u = eigenvec_t1_u(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_t1_u = atan2(largest_eigenvec_t1_u(2), largest_eigenvec_t1_u(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_t1_u < 0)
    angle_t1_u = angle_t1_u + 2*pi;
end

phi_t1_u = angle_t1_u;
a_t1_u = chisquare_val * sqrt(largest_eigenval_t1_u);
b_t1_u = chisquare_val * sqrt(smallest_eigenval_t1_u);

% the ellipse in x and y coordinates 
ellipse_x_r_t1_u  = a_t1_u*cos(theta_grid);
ellipse_y_r_t1_u  = b_t1_u*sin(theta_grid);

% rotation matrix of angle phi
rot_t1_u = [cos(phi_t1_u) sin(phi_t1_u); -sin(phi_t1_u) cos(phi_t1_u)];

% rotation of the ellipse to angle phi
r_ellipse_t1_u = [ellipse_x_r_t1_u;ellipse_y_r_t1_u]' * rot_t1_u;

%% figure 1

fig1 = figure(1);
ax1 = axes('Parent',fig1);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
axis equal

plt_t0(1) = plot(x_t0,y_t0,'Marker','.','LineStyle','none','Color','#0072BD');
plt_t0(2) = plot(r_ellipse_t0(:,1) + mu_t0(1),r_ellipse_t0(:,2) + mu_t0(2),'Color','#0072BD');
plt_t0(3) = plot(mu_t0(1),mu_t0(2),'Marker','*','LineStyle','none','Color','#0072BD','LineWidth',2);

plt_t1_pst(1) = plot(x_vec_t1(1,:),x_vec_t1(2,:),'Marker','.','LineStyle','none','Color','#77AC30');
plt_t1_pst(2) = plot(r_ellipse_t1_psb(:,1) + mu_t1_mc(1),r_ellipse_t1_psb(:,2) + mu_t1_mc(2),'Color','#77AC30');
plt_t1_pst(3) = plot(mu_t1_mc(1),mu_t1_mc(2),'Marker','*','LineStyle','none','Color','#77AC30','LineWidth',2);

% h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 2

fig2 = figure(2);
ax2 = axes('Parent',fig2);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
% axis equal

copyobj(plt_t0(2:3),ax2)

plt_t1_bar(2) = plot(r_ellipse_t1_bar(:,1) + mu_t1_bar(1),r_ellipse_t1_bar(:,2) + mu_t1_bar(2),'Color','#A2142F');
plt_t1_bar(3) = plot(mu_t1_bar(1),mu_t1_bar(2),'Marker','*','LineStyle','none','Color','#A2142F','LineWidth',2);

% h = legend([plt_t1_bar,plt_t1_act],'$\overline{bel}(x_1)$','$\overline{UE}_{1,95\%}$','$\bar{\mu}_1$','$x_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 3

fig3 = figure(3);
ax3 = axes('Parent',fig3);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
axis equal

copyobj(plt_t0(2),ax3)
plt_t0(4) = plot(chi_t0(1,[1:3 5:6]),chi_t0(2,[1:3 5:6]),'Marker','*','LineStyle','none','Color','#0072BD','LineWidth',2);

plt_t1_bar_u(1) = plot(r_ellipse_t1_bar_u(:,1) + mu_t1_bar_u(1),r_ellipse_t1_bar_u(:,2) + mu_t1_bar_u(2),'Color','#7E2F8E');
plt_t1_bar_u(2) = plot(chi_t1_bar(1,[1:3 5:6]),chi_t1_bar(2,[1:3 5:6]),'Marker','*','LineStyle','none','Color','#7E2F8E','LineWidth',2);

% h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 4

fig4 = figure(4);
ax4 = axes('Parent',fig4);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
axis equal

copyobj(plt_t1_pst(2:3),ax4)
copyobj(plt_t1_bar(2:3),ax4)
copyobj(plt_t1_bar_u(1),ax4)
plt_t1_bar_u(3) = plot(mu_t1_bar_u(1),mu_t1_bar_u(2),'Marker','*','LineStyle','none','Color','#7E2F8E','LineWidth',2);

% h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 5

fig5 = figure(5);
ax5 = axes('Parent',fig5);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
axis equal

copyobj(plt_t1_pst(2:3),ax5)

plt_t1_act = plot(x_t1(1),x_t1(2),'*k','LineWidth',2);

% plt_z_t1_vec = plot(z_t1_vec,zeros(size(z_t1_vec)),'Marker','.','LineStyle','none','Color','#00FFFF');

plt_z_t1(1) = xline(z_t1,'LineStyle','--','Color','#EDB120');
plt_z_t1(2) = xline(z_t1-chisquare_val*sqrt(Q),'Color','#EDB120');
plt_z_t1(3) = xline(z_t1+chisquare_val*sqrt(Q),'Color','#EDB120');

% plot the error ellipse
plt_t1_kf(1) = plot(r_ellipse_t1_kf(:,1) + mu_t1_kf(1),r_ellipse_t1_kf(:,2) + mu_t1_kf(2),'-g');
plt_t1_kf(2) = plot(mu_t1_kf(1),mu_t1_kf(2),'*g','LineWidth',2);

% plt_t1_pst(1) = plot(x_vec_t1_pst(1,:),x_vec_t1_pst(2,:),'.g');
% plt_t1_pst(2) = plot(r_ellipse_t1_psb(:,1) + mu_t1_psb(1),r_ellipse_t1_psb(:,2) + mu_t1_psb(2),'-g');
% plt_t1_pst(3) = plot(mu_t1_psb(1),mu_t1_psb(2),'*g','LineWidth',2);

% h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 6

fig6 = figure(6);
ax6 = axes('Parent',fig6);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
% axis equal

copyobj(plt_t1_bar(2:3),ax6)
copyobj(plt_t1_act,ax6)
copyobj(plt_z_t1,ax6)

% plot the error ellipse
plt_t1(1) = plot(r_ellipse_t1(:,1) + mu_t1(1),r_ellipse_t1(:,2) + mu_t1(2),'-r');
plt_t1(2) = plot(mu_t1(1),mu_t1(2),'*r','LineWidth',2);

% h = legend([plt_t1_bar,plt_t1_act],'$\overline{bel}(x_1)$','$\overline{UE}_{1,95\%}$','$\bar{\mu}_1$','$x_1$');
% set(h,'interpreter','Latex','FontSize',12)

%% figure 7

fig7 = figure(7);
ax7 = axes('Parent',fig7);
grid on
hold on
xlabel('x')
ylabel('y')
xlim([-2, 3]);
% ylim([-2, 2]);
axis equal

copyobj(plt_t1_bar_u([1 3]),ax7)

copyobj(plt_t1_act,ax7)
copyobj(plt_z_t1,ax7)

% plot the error ellipse
plt_t1_u(1) = plot(r_ellipse_t1_u(:,1) + mu_t1_u(1),r_ellipse_t1_u(:,2) + mu_t1_u(2),'-m');
plt_t1_u(2) = plot(mu_t1_u(1),mu_t1_u(2),'*m','LineWidth',2);

% h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% set(h,'interpreter','Latex','FontSize',12)

% %% figure 4
% 
% fig4 = figure(4);
% ax4 = axes('Parent',fig4);
% grid on
% hold on
% xlabel('x')
% ylabel('y')
% xlim([-2, 3]);
% % ylim([-2, 2]);
% axis equal
% 
% copyobj(plt_t1_pst(2:3),ax4)
% copyobj(plt_t1_bar(2:3),ax4)
% copyobj(plt_t1_u(1),ax4)
% plot(mu_t1_bar_u(1),mu_t1_bar_u(2),'*m','LineWidth',2);
% 
% % h = legend('$bel(x_0)$','$UE_{0,95\%}$','$\mu_0$','$\overline{bel}(x`_1)$','$\overline{UE}`_{1,95\%}$','$\bar{\mu}`_1$');
% % set(h,'interpreter','Latex','FontSize',12)
