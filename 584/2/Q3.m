% Q4
% 2D Case 2L and 1Star
% 1 subtended angle and 1 bearing
clc
clear all
close all
format long

rad = pi/180;

% P0: my location
% psi L2 to star
% theta subtended angle


L1 = [0 0 0];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
z1 = L1(3);

L2 = [5 0 0];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
z2 = L2(3);

theta = 90*rad;      
psi   = 135*rad;    

x = [x1 x2];
y = [y1 y2];
z = [z1 z2];
XS   = [0 1 0]';
% crate a grid of initial points
initXLim    = [-5:2:5];   % range for initial point x values
initYLim    = [-5:2:5];   % range for initial point y values
initZLim    = [-5:2:5];   % range for initial point z values

init_matrix = [];

for ai=1:length(initXLim)
    for aj=1:length(initYLim)
        for ak=1:length(initZLim)
        init = [initXLim(ai) initYLim(aj) initZLim(ak)];
        init_matrix = [init_matrix; init];
        end
    end
end

% optimization
F       = @(x) HW2Q2costFun(x,L1,L2,theta,psi,XS);                   % define the cost function
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12,'Display','none');    % set the tolerance

for ck=1:length(init_matrix)
    [xmin(ck,:) , fval(ck,:)] = fminunc(F, init_matrix(ck,:), options);
end

% find the distance between fixes and lighthouses
fixL1Dist = sqrt((xmin(:,1)-x(1)).^2 + (xmin(:,2)-y(1)).^2 + (xmin(:,3)-z(1)).^2) ; % distance to L1 at each fix
fixL2Dist = sqrt((xmin(:,1)-x(2)).^2 + (xmin(:,2)-y(2)).^2 + (xmin(:,3)-z(2)).^2) ; % distance to L2 at each fix
xmin( find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),: ) = [];

x = xmin(:,1);
y = xmin(:,2);
z = xmin(:,3);

P0 = [2.5 2.5 0];
%% plots
%2D projections

figure(5)

subplot(3,1,1)
scatter(P0(1), P0(3),50,'r','fill')
hold on
scatter(x,z,50,'b','fill')
hold off
grid on
ylabel('z [m]')
xlabel('x [m]')
legend('Actual Location (xz)','Position Fixes','Location','best')
title('Problem 3: 2D Plot XZ','2 Lighthouses and 1 Star')


subplot(3,1,2)
scatter(P0(1), P0(2),50,'r','fill')
hold on
scatter(x,y,50,'b','fill')
scatter(P0(1), P0(2),50,'r','fill')
hold off
grid on
ylabel('y [m]')
xlabel('x [m]')
legend('Actual Location (xy)','Position Fixes','Location','best')
title('Problem 3: 2D Plot XY','2 Lighthouses and 1 Star')

subplot(3,1,3)
scatter(z,y,50,'b','fill')
hold on
scatter(P0(3), P0(2),50,'r','fill')
hold off
ylabel('y [m]')
xlabel('z [m]')
grid on
legend('Actual Location (zy)','Position Fixes','Location','best')
title('Problem 3: 2D Plot ZY','2 Lighthouses and 1 Star')

figure(6)
Xsp=2.5;
Ysp=2.5;
Zsp=2.5;
[Xsp,Ysp,Zsp] = sphere;
surf(2.5.*Xsp + 2.5,2.5.*Ysp,2.5.*Zsp,'FaceAlpha',0.25);
hold on
scatter3(x,y,z,50,'b','fill')
scatter3(P0(1), P0(2), P0(3),50,'r','fill')
hold off
zlabel('z [m]')
ylabel('y [m]')
xlabel('x [m]')
title('Problem 3: 3D Plot XYZ','2 Lighthouses and 1 Star')

