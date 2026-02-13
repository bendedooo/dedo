% Q4
% 3D Case 2L and 2Star
% 1 subtended angle and 2 bearing
clc
clear all
close all
format long
warning('off')
rad = pi/180;

% P0: my location
% psi to star
% theta subtended angle


L1 = [0 0 0];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
z1 = L1(3);

L2 = [5 0 0];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
z2 = L2(3);

theta = 90*rad;     % measurement from L2
psi1  = 135*rad;   % measurement from S1 to L2
psi2  = 90*rad;    % measurement from S2 to L2

x = [x1 x2];
y = [y1 y2];
z = [z1 z2];

XS1   = [0 1 0]';
XS2   = [0 0 1]';

% crate a grid of initial points
initXLim    = [-5:3:5];   % range for initial point x values
initYLim    = [-5:3:5];   % range for initial point y values
initZLim    = [-5:3:5];   % range for initial point z values

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
F       = @(x) HW2Q4costFun(x,L1,L2,theta,psi1,psi2,XS1,XS2);                   % define the cost function
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12,'Display','none');    % set the tolerance

for ck=1:length(init_matrix)
    [xmin1(ck,:) , fval1(ck,:)] = fminunc(F, init_matrix(ck,:), options);
end

% find the distance between fixes and lighthouses
fixL1Dist = ((xmin1(:,1)-x(1)).^2 + (xmin1(:,2)-y(1)).^2 + (xmin1(:,3)-z(1)).^2).^0.5 ; % distance to L1 at each fix
fixL2Dist = ((xmin1(:,1)-x(2)).^2 + (xmin1(:,2)-y(2)).^2 + (xmin1(:,3)-z(2)).^2).^0.5 ; % distance to L2 at each fix

xmin1( find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),: ) = [];
fval1( find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),: ) = [];


ind1 = find(fval1==min(fval1));

xout = xmin1(ind1,1);
yout = xmin1(ind1,2);
zout = xmin1(ind1,3);

P0 = [2.5 2.5 0];


%% plots
%2D projections

figure(7)

subplot(3,1,1)
scatter(x,z,50,'g','fill')
hold on
scatter(xout,zout,50,'r','fill')
hold off
grid on
ylabel('z [m]')
xlabel('x [m]')
legend('Lighthouse Locations (xz)','Position Fix','Location','best')
title('Problem 4: 2D Plot XZ','3D Case Two Lighthouses and Two Stars')


subplot(3,1,2)
scatter(x,y,50,'g','fill')
hold on
scatter(xout,yout,50,'r','fill')
hold off
grid on
ylabel('y [m]')
xlabel('x [m]')
legend('Lighthouse Locations (xy)','Position Fix','Location','best')
title('Problem 4: 2D Plot XY','3D Case Two Lighthouses and Two Stars')

subplot(3,1,3)
scatter(z,y,50,'g','fill')
hold on
scatter(zout,yout,50,'r','fill')
hold off
ylabel('y [m]')
xlabel('z [m]')
grid on
legend('Lighthouse Locations (zy)','Position Fix','Location','best')
title('Problem 4: 2D Plot ZY','3D Case Two Lighthouses and Two Stars')

figure(8)
Xsp=2.5;
Ysp=2.5;
Zsp=2.5;
[Xsp,Ysp,Zsp] = sphere;
surf(2.5.*Xsp + 2.5,2.5.*Ysp,2.5.*Zsp,'FaceAlpha',0.25);
hold on
scatter3(P0(1), P0(2), P0(3),50,'b','fill')
scatter3(xout, yout, zout,50,'r','fill')
scatter3(x,y,z,50,'g','fill')
hold off
zlabel('z [m]')
ylabel('y [m]')
xlabel('x [m]')
legend('','Lighthouse Locations','Position Fix','Actual Location','Location','best')
title('Problem 4: 3D Plot XYZ','3D Case Two Lighthouses and Two Stars')

