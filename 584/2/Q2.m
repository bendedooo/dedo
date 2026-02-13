clc
clear all
close all
format long

rad = pi/180;

%% PART A
% 2D Position Fixing with 1 Subtended Angle and 1 Bearing

L1 = [0 0];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
theta1 = -165*rad;    %  measurement from L1

L2 = [4 2];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
theta2 = 150*rad;      % measurement from L2

theta = 45*rad;
psi   = 150*rad;

x  = [x1 x2];
y  = [y1 y2];
XS = [0 1]';

% crate a grid of initial points
initXLim    = [-2:2:10];   % range for initial point x values
initYLim    = [-2:2:10];   % range for initial point y values
init_matrix = [];          % matrix of initial points

for ai=1:length(initXLim)
    for aj=1:length(initYLim)
        init = [initXLim(ai) initYLim(aj)];
        init_matrix = [init_matrix; init]; 
    end
end

F       = @(x) HW2Q2costFun(x,L1,L2,theta,psi,XS);                   % define the cost function

% optimization
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12,'Display','none');    % set the tolerance

 
for ck=1:length(init_matrix)
    [xmin(ck,:) , fval(ck,:)] = fminunc(F, init_matrix(ck,:), options);
end



% find the distance between fixes and lighthouses
fixL1Dist = ((xmin(:,1)-x(1)).^2 + (xmin(:,2)-y(1)).^2).^(1/2) ; % distance to L1 at each fix
fixL2Dist = ((xmin(:,1)-x(2)).^2 + (xmin(:,2)-y(2)).^2).^(1/2) ; % distance to L2 at each fix
xmin(find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),:) = [];
fval(find(1-(fixL1Dist>0.01).*(fixL2Dist>0.01)),:) = [];
ind = find(fval==min(fval));

figure(3)
scatter(x,y,40,'g','fill')
hold on
scatter(xmin(ind,1),xmin(ind,2),40,'r','fill')
% scatter(xmin(:,1),xmin(:,2),10,'m','fill')
title('Problem 2: Part a','2D Position Fix with 1 Subtended Angle and 1 Bearing Measurement')
xlabel('x [m]')
ylabel('y [m]')
legend('Lighthouses','Position Fix','Location','best')
grid on

%% PART B

L1 = [0 0];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
theta1 = -140*rad;    

L2 = [4 2];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
theta2 = 90*rad;      

L3 = [1 4];  
x3 = L3(1);
y3 = L3(2);
theta3 = -30*rad;      


x = [x1 x2 x3];
y = [y1 y2 y3];

% L1 and L2
psi12   = theta2;
theta12 = 2*pi - (-theta1+theta2); %130*rad
F       = @(x) HW2Q2costFun(x,L1,L2,theta12,psi12,XS);                   % define the cost function
for i12=1:length(init_matrix)
    [xmin12(i12,:) , fval12(i12,:)] = fminunc(F, init_matrix(i12,:), options);
end


% L1 and L3 
psi13   = theta3;
theta13 = (theta3-theta1);
F       = @(x) HW2Q2costFun(x,L1,L3,theta13,psi13,XS);                   % define the cost function
for i13=1:length(init_matrix)
    [xmin13(i13,:) , fval13(i13,:)] = fminunc(F, init_matrix(i13,:), options);
end


% L2 and L3
psi23   = theta2;
theta23 = (theta2-theta3);
F       = @(x) HW2Q2costFun(x,L3,L2,theta23,psi23,XS);                   % define the cost function
for i23=1:length(init_matrix)
    [xmin23(i23,:) , fval23(i23,:)] = fminunc(F, init_matrix(i23,:), options);
end




% find the distance between fixes and lighthouses
fix12L1Dist = sqrt((xmin12(:,1)-x(1)).^2 + (xmin12(:,2)-y(1)).^2) ; % distance to L1 at each fix
fix12L2Dist = sqrt((xmin12(:,1)-x(2)).^2 + (xmin12(:,2)-y(2)).^2) ; % distance to L2 at each fix
xmin12( find(1-(fix12L1Dist>0.01).*(fix12L2Dist>0.01)),: ) = [];
fval12( find(1-(fix12L1Dist>0.01).*(fix12L2Dist>0.01)),: ) = [];

ind12 = find(fval12==min(fval12));

fix13L1Dist = sqrt((xmin13(:,1)-x(1)).^2 + (xmin13(:,2)-y(1)).^2) ; % distance to L1 at each fix
fix13L2Dist = sqrt((xmin13(:,1)-x(2)).^2 + (xmin13(:,2)-y(2)).^2) ; % distance to L2 at each fix
xmin13( find(1-(fix13L1Dist>0.01).*(fix13L2Dist>0.01)),: ) = [];
fval13( find(1-(fix13L1Dist>0.01).*(fix13L2Dist>0.01)),: ) = [];

ind13 = find(fval13==min(fval13));

fix23L1Dist = sqrt((xmin23(:,1)-x(1)).^2 + (xmin23(:,2)-y(1)).^2) ; % distance to L1 at each fix
fix23L2Dist = sqrt((xmin23(:,1)-x(2)).^2 + (xmin23(:,2)-y(2)).^2) ; % distance to L2 at each fix
xmin23( find(1-(fix23L1Dist>0.01).*(fix23L2Dist>0.01)),: ) = [];
fval23( find(1-(fix23L1Dist>0.01).*(fix23L2Dist>0.01)),: ) = [];

ind23 = find(fval23==min(fval23));
xmins = [xmin12(ind12,1), xmin13(ind13,1), xmin23(ind23,1) ; xmin12(ind12,2), xmin13(ind13,2), xmin23(ind23,2) ];

figure(4)
scatter(x,y,40,'g','fill')
hold on
% scatter(xmin12(ind,1),xmin12(ind,2),40,'b','fill')
% scatter(xmin13(ind,1),xmin13(ind,2),40,'b','fill')
% scatter(xmin23(ind,1),xmin23(ind,2),40,'b','fill')
% scatter(xmin12(:,1),xmin12(:,2),40,'r','fill')
% scatter(xmin13(:,1),xmin13(:,2),40,'r','fill')
% scatter(xmin23(:,1),xmin23(:,2),40,'r','fill')
scatter(xmins(1,:),xmins(2,:),'b','fill')
scatter(sum(xmins(1,:))/length(xmins(1,:)),sum(xmins(2,:))/length(xmins(2,:)),40,'r','fill')
title('Problem 2: Part b','3 Lighthouses')
xlabel('x [m]')
ylabel('y [m]')
legend('Lighthouses','Position Fixes','Lighthouse Center','Location','best')
grid on


