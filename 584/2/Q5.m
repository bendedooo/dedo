% Q5
% 3D Case 2L and 2Star noise added
% 1 subtended angle and 2 bearing
clc
clear all
close all
format long
warning('off')
rad = pi/180;

L1 = [0 0 0];  % first known location coordinates [x y z]
x1 = L1(1);
y1 = L1(2);
z1 = L1(3);

L2 = [5 0 0];  % second known location coordinates [x y z]
x2 = L2(1);
y2 = L2(2);
z2 = L2(3);

theta = 90*rad;    % measurement from L2
psi1  = 135*rad;   % measurement from S1 to L2
psi2  = 90*rad;    % measurement from S2 to L2

seed = 1333;
s = rng;
s.Seed = seed;
rng(s.Seed);

runNumber = 200;
measNoise = (2*randn(3,200))*rad;
thetaMatrix = theta + measNoise(1,:);
psi1Matrix = psi1 + measNoise(2,:);
psi2Matrix = psi2 + measNoise(3,:);

x = [x1 x2];
y = [y1 y2];
z = [z1 z2];

XS1   = [0 1 0]';
XS2   = [0 0 1]';

% crate a grid of initial points
initXLim    = [-1:3:5];   % range for initial point x values
initYLim    = [-1:3:5];   % range for initial point y values
initZLim    = [-1:3:5];   % range for initial point z values

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

fixMatrix = [];

for i=1:runNumber
    F       = @(x) HW2Q4costFun(x,L1,L2,thetaMatrix(1,i),psi1Matrix(1,i),psi2Matrix(1,i),XS1,XS2);                   % define the cost function
    options = optimoptions(@fminunc,'OptimalityTolerance',1e-12,'Display','none');    % set the tolerance

    for ck=1:length(init_matrix)
        [xmin1(ck,:) , fval1(ck,:)] = fminunc(F, init_matrix(ck,:), options);
    end

    % find the distance between fixes and lighthouses
    fixL1Dist = ((xmin1(:,1)-x(1)).^2 + (xmin1(:,2)-y(1)).^2 + (xmin1(:,3)-z(1)).^2).^0.5 ; % distance to L1 at each fix
    fixL2Dist = ((xmin1(:,1)-x(2)).^2 + (xmin1(:,2)-y(2)).^2 + (xmin1(:,3)-z(2)).^2).^0.5 ; % distance to L2 at each fix

    xmin1( find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),: ) = [];
    fval1( find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),: ) = [];


    % take the min cost

    ind1 = find(fval1==min(fval1));
    xout = xmin1(ind1,1);
    yout = xmin1(ind1,2);
    zout = xmin1(ind1,3);

    fixMatrix = [fixMatrix; xout yout zout];

end

P0 = [2.5 2.5 0];

coordError = fixMatrix-((ones(length(fixMatrix),1)).*P0);
posError = vecnorm(coordError');

%% plots
%2D projections

figure(9)

subplot(3,1,1)
scatter(x,z,50,'g','fill')
hold on
scatter(fixMatrix(:,1),fixMatrix(:,3),50,'b','fill')
scatter(P0(1),P0(3),50,'r','fill')
hold off
grid on
ylabel('z [m]')
xlabel('x [m]')
legend('Lighthouse Locations (xz)','Position Fixes','Actual Location','Location','best')
title('Problem 5: 2D Plot XZ','3D Case Two Lighthouses and Two Stars Noisy Measurements')


subplot(3,1,2)
scatter(x,y,50,'g','fill')
hold on
scatter(fixMatrix(:,1),fixMatrix(:,2),50,'b','fill')
scatter(P0(1),P0(2),50,'r','fill')
hold off
grid on
ylabel('y [m]')
xlabel('x [m]')
legend('Lighthouse Locations (xy)','Position Fixes','Actual Location','Location','best')
title('Problem 5: 2D Plot XY','3D Case Two Lighthouses and Two Stars Noisy Measurements')

subplot(3,1,3)
scatter(z,y,50,'g','fill')
hold on
scatter(fixMatrix(:,3),fixMatrix(:,2),50,'b','fill')
scatter(P0(3),P0(2),50,'r','fill')
hold off
ylabel('y [m]')
xlabel('z [m]')
grid on
legend('Lighthouse Locations (zy)','Position Fixes','Actual Location','Location','best')
title('Problem 5: 2D Plot ZY','3D Case Two Lighthouses and Two Stars Noisy Measurements')

figure(10)
Xsp=2.5;
Ysp=2.5;
Zsp=2.5;
[Xsp,Ysp,Zsp] = sphere;
surf(2.5.*Xsp + 2.5,2.5.*Ysp,2.5.*Zsp,'FaceAlpha',0.25)
hold on
scatter3(x,y,z,50,'g','fill')
scatter3(fixMatrix(:,1), fixMatrix(:,2),fixMatrix(:,3),50,'b','fill')
scatter3(P0(1), P0(2), P0(3),50,'r','fill');
hold off
zlabel('z [m]')
ylabel('y [m]')
xlabel('x [m]')
legend('','Lighthouse Locations','Position Fixes','Actual Location','Location','best')
title('Problem 5: 3D Plot XYZ','3D Case Two Lighthouses and Two Stars Noisy Measurements')

figure(11)
hist(posError,20)
xlabel('Position Error [m]')
ylabel('Samples')
title('Problem 5: Position Errors','3D Case Two Lighthouses and Two Stars Noisy Measurements')
