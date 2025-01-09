clc, clear all, close all

% robot pose
xR_I = 0;
yR_I = 0;
thetaR_I = 13.5*pi/4;

thetaR_I = normalizeFunc(thetaR_I,-pi,pi);
    
R_pose = [xR_I;yR_I;thetaR_I];

% sensor parameters
rng = 0.5;
psi = pi/1.2;
res = pi/20;

S_para = [rng;psi;res]; 

% sensor pose
dS_R = [0.05;0]; % sensor distance from robot center
thetaS_R = 0;

xS_I = xR_I + dS_R(1)*cos(thetaR_I) - dS_R(2)*sin(thetaR_I);
yS_I = yR_I + dS_R(1)*sin(thetaR_I) + dS_R(2)*cos(thetaR_I);
thetaS_I = thetaR_I + thetaS_R;

S_pose = [xS_I;yS_I;thetaS_I];
    
% sesor error parameters
sigma_d = 0;
sigma_phi = 0;
sigma_s = 0;

% landmark positions
% marker(1) = 1;
% point(1).x = 0.18;
% point(1).y = -0.2;
% point(1).s = '^';

marker(1) = 1;
point(1).x = -0.2;
point(1).y = -0.275;
point(1).s = 'o';

marker(2) = 1;
point(2).x = -0.2;
point(2).y = 0.175;
point(2).s = 's';

marker(3) = 1;
point(3).x = 0.25;
point(3).y = -0.275;
point(3).s = '^';

marker(4) = 1;
point(4).x = 0.28-0.1;
point(4).y = -0.29+0.1;
point(4).s = 'p';

% marker(3) = 1;
% point(3).x = 0.2;
% point(3).y = 0;
% point(3).s = '^';
% 
% marker(4) = 1;
% point(4).x = -0.16;
% point(4).y = -0.2;
% point(4).s = '^';

marker(5) = 1;
point(5).x = 0.4;
point(5).y = 0.07;
point(5).s = 'o';

% marker(5) = 1;
% point(5).x = 0.18;
% point(5).y = -0.2;
% point(5).s = '^';

nP = sum(marker);
% nP = length(marker);

xM = zeros(1,nP);
yM = zeros(1,nP);
sM = zeros(1,nP);
k = 0;
for i = 1:length(marker)
    if marker(i) == 1
        k = k+1;
        xM(k) = point(i).x;
        yM(k) = point(i).y;
        sM(k) = point(i).s;
    end
end

% noise distribution type
app.dist_type = 0;

% sensor model
rti = zeros(1,nP);
phi = zeros(1,nP);
sti = zeros(1,nP);
for j = 1:nP
    rti(j) = sqrt((xM(j)-xS_I)^2 + (yM(j)-yS_I)^2) + sample(app,sigma_d^2);

    phi(j) = atan2(yM(j)-yS_I, xM(j)-xS_I) - thetaR_I + sample(app,sigma_phi^2);    
    phi(j) = normalizeFunc(phi(j),-pi,pi);

    sti(j) = markerS2I(char(sM(j)));    
    sti(j) = round(sti(j) + sample(app,sigma_s^2));
end

F_temp = [rti;phi;sti];

ang = floor(psi/res);

ind = zeros(1,nP);
phiD = zeros(1,nP);
for i = 1:1:2*ang+1
    for j = 1:nP
        if rti(1,j) < rng && phi(j) > res*(i-1-ang) && phi(j) < res*(i-ang)
            if j > 1
                if rti(1,j-1) < rng && phi(j-1) > res*(i-1-ang) && phi(j-1) < res*(i-ang)
                    if rti(1,j) < rti(1,j-1)
                        ind(j-1) = 0;
                        ind(j) = i;
                    end
                else
                    ind(j) = i;
                end
            else
                ind(j) = i;
            end
        end
        phiD(j) = thetaS_I+res*(ind(j)-ang-0.5);
    end
end

F_info = [rti;phiD;ind];

d_i = rng*ones(1,2*ang);
mrk = zeros(1,2*ang);

for j = 1:length(rti)
    if ind(j) > 0
        d_i(ind(j)) = rti(j);
        mrk(ind(j)) = sti(j);
    end
end

mx = zeros(1,nP);
my = zeros(1,nP);
for i = 1:nP
    mx(1,i) = xS_I+rti(i)*cos(thetaS_I+res*(ind(i)-ang-0.5));
    my(1,i) = yS_I+rti(i)*sin(thetaS_I+res*(ind(i)-ang-0.5));
end

M = [mx;my;sti];

figure
hold(gca,'on')
grid(gca,'on')
axis(gca,'equal')
axis(gca,'padded')

drawRobot(gca,R_pose,S_pose,S_para,F_info,M)

for i = 1:nP
    plot(gca,xM(i),yM(i),char(sM(i)),'MarkerEdgeColor','k')
end

function drawRobot(ax,R_pose,S_pose,S_para,F_info,M)

    rR = 0.1;
    rS = 0.025;

    xR = R_pose(1);
    yR = R_pose(2);
    tR = R_pose(3);

    hold(ax,"on")
    plot(ax,xR + rR*cos(linspace(0,2*pi)),yR + rR*sin(linspace(0,2*pi)),'-k',xR,yR,'.k')
    plot(ax,[xR,xR + rR*cos(tR)],[yR,yR + rR*sin(tR)],'-k')

    rng = S_para(1);
    psi = S_para(2);
    res = S_para(3);

    xS = S_pose(1);
    yS = S_pose(2);
    tS = S_pose(3);

    ang = floor(psi/res);
    plot(ax,xS + rS*cos(linspace(0,2*pi)),yS + rS*sin(linspace(0,2*pi)),'-b',xS,yS,'.r')
    pcd = plot(ax,xS + rng*cos(linspace(tS-psi,tS+psi)),yS + rng*sin(linspace(tS-psi,tS+psi)),':r');
    pcd.Color = [pcd.Color 0.25];
    pnd = plot(ax,[xS xS + rng*cos(tS-psi)],[yS yS + rng*sin(tS-psi)],':r');
    pnd.Color = [pnd.Color 0.25];
    ppd = plot(ax,[xS xS + rng*cos(tS+psi)],[yS yS + rng*sin(tS+psi)],':r');
    ppd.Color = [ppd.Color 0.25];
    for i = 1:2*ang+1
        pmd = plot(ax,[xS, xS+rng*cos(tS+res*(i-1-ang))],[yS, yS+rng*sin(tS+res*(i-1-ang))],':r');
        pmd.Color = [pmd.Color 0.25];
    end
    
    dF = F_info(1,:);
    tF = F_info(2,:);
    iF = F_info(3,:);
%     sF = F_info(3,:);
    
    d_i = rng*ones(1,2*ang);
    t_i = zeros(1,2*ang);
    mx = zeros(1,2*ang);
    my = zeros(1,2*ang);
    ms = zeros(1,2*ang);

    for j = 1:length(dF)
        if iF(j) > 0
            d_i(iF(j)) = dF(j);
            t_i(iF(j)) = tF(j);
            mx(iF(j)) = M(1,j);
            my(iF(j)) = M(2,j);
            ms(iF(j)) = M(3,j);
        end
    end
    
    for i = 1:2*ang
        if d_i(i) < rng
%             plot(ax, xS+d_i(i)*cos(linspace(tS+res*(i-1-ang),tS+res*(i-ang))), yS + d_i(i)*sin(linspace(tS+res*(i-1-ang),tS+res*(i-ang))), '-r')
%             plot(ax, [xS, xS+d_i(i)*cos(tS+res*(i-1-ang))], [yS, yS+d_i(i)*sin(tS+res*(i-1-ang))], '-r')
%             plot(ax, [xS, xS+d_i(i)*cos(tS+res*(i-0.5-ang))], [yS, yS+d_i(i)*sin(tS+res*(i-0.5-ang))],'--r')
%             plot(ax, [xS, xS+d_i(i)*cos(tS+res*(i-ang))], [yS, yS+d_i(i)*sin(tS+res*(i-ang))], '-r')
%             plot(ax, xS+d_i(i)*cos(tS-ang*res+(i*res)-res/2), yS+d_i(i)*sin(tS-ang*res+(i*res)-res/2),'r','Marker',markerI2S(mrk(i)))
            plot(ax, xS+d_i(i)*cos(linspace(t_i(i)-0.5*res,t_i(i)+0.5*res)), yS+d_i(i)*sin(linspace(t_i(i)-0.5*res,t_i(i)+0.5*res)), '-r')
            plot(ax, [xS, xS+d_i(i)*cos(t_i(i)-0.5*res)], [yS, yS+d_i(i)*sin(t_i(i)-0.5*res)], '-r')
            plot(ax, [xS, xS+d_i(i)*cos(t_i(i))], [yS, yS+d_i(i)*sin(t_i(i))],'--r')
            plot(ax, [xS, xS+d_i(i)*cos(t_i(i)+0.5*res)], [yS, yS+d_i(i)*sin(t_i(i)+0.5*res)], '-r')
            plot(ax, mx(i), my(i), 'r', 'Marker', markerI2S(ms(i)))
        else
            pn = plot(ax, [xS, xS+rng*cos(tS+res*(i-1-ang))], [yS, yS+rng*sin(tS+res*(i-1-ang))], '-r');
            pn.Color = [pn.Color 0.25];
            pp = plot(ax, [xS, xS+rng*cos(tS+res*(i-ang))], [yS, yS+rng*sin(tS+res*(i-ang))], '-r');
            pp.Color = [pp.Color 0.25];
        end
    end
    
end

function n = normalizeFunc(value,start,finish) 

  width       = finish - start   ;
  offsetValue = value - start ;

  n = ( offsetValue - ( floor( offsetValue / width ) * width ) ) + start ;
  
end

function mS = markerI2S(mI)
    if rem(mI,4) == 1
        mS = 's';
    elseif rem(mI,4) == 2
        mS = '^';
    elseif rem(mI,4) == 3
        mS = 'p';
    else
        mS = 'o';
    end
end

function mI = markerS2I(mS)
    if mS == 's'
        mI = 1;
    elseif mS == '^'
        mI = 2;
    elseif mS == 'p'
        mI = 3;
    else
        mI = 0;
    end
end
