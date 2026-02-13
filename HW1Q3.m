% PROBLEM 3
clc
clear all
close all
format long

% beacons and measurements
L1 = [2 2];    % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
R1 = 4.8;      % range measurement from L1

L2 = [6 6];    % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
R2 = 8.4;        % range measurement from L2

L3 = [10 5];  % third known location coordinates [x y]
x3 = L3(1);
y3 = L3(2);
R3 = 9.21;        % range measurement from L3

% initial position guess
init1 = [14    2]; 
init2 = [-2   16];
init3 = [0.6 4.8];

F = @(x) HW1Q3costFun(x,L1,L2,L3,R1,R2,R3);                   % define the cost function
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12); % set the tolerance

%% PART A

%optimization
[xmin1, fval1] = fminunc(F, init1, options);
[xmin2, fval2] = fminunc(F, init2, options);
[xmin3, fval3] = fminunc(F, init3, options);

figure(1)
viscircles(L1,R1,'Color','g');
viscircles(L2,R2,'Color','r');
viscircles(L3,R3,'Color','b');
hold on
scatter(xmin1(1),xmin1(2),40,'m','fill')
scatter(xmin2(1),xmin2(2),40,'c','fill')
scatter(xmin3(1),xmin3(2),40,'k','fill')
scatter(init1(1),init1(2),40,'m*')
scatter(init2(1),init2(2),40,'c*')
scatter(init3(1),init3(2),40,'k*')
title('Problem 3: Part a','2D Position Fix with Three Range Measurements')
xlabel('x')
ylabel('y')
legend('Fix 1','Fix 2','Fix 3','Initial Guess 1','Initial Guess 2','Initial Guess 3','Location','best')
grid on

save ('Q3','xmin1','L1','L2','L3','R1','R2','R3')

%% PART B 

% measurements for inconsistent case
R1 = 7;
R2 = 8;
R3 = 6;

% initial position guess
init1 = [12    3]; 
init2 = [-5   18];
init3 = [0.2 8.8];

% optimization
F = @(x) HW1Q3costFun(x,L1,L2,L3,R1,R2,R3);  % new cost function

[xmin1, fval1] = fminunc(F, init1, options);
[xmin2, fval2] = fminunc(F, init2, options);
[xmin3, fval3] = fminunc(F, init3, options);

figure(2)
viscircles(L1,R1,'Color','g');
viscircles(L2,R2,'Color','r');
viscircles(L3,R3,'Color','b');
hold on
scatter(xmin1(1),xmin1(2),40,'m','fill')
scatter(xmin2(1),xmin2(2),40,'c','fill')
scatter(xmin3(1),xmin3(2),40,'k','fill')
scatter(init1(1),init1(2),40,'m*')
scatter(init2(1),init2(2),40,'c*')
scatter(init3(1),init3(2),40,'k*')
title('Problem 3: Part b','Inconsistent 2D Position Fix with Three Range Measurements')
xlabel('x')
ylabel('y')
legend('Fix 1','Fix 2','Fix 3','Initial Guess 1','Initial Guess 2','Initial Guess 3','Location','best')
grid on

%% PART C

% crate a grid of initial points
initXLim    = [-8:2:24];   % range for initial point x values
initYLim    = [-8:3:24];   % range for initial point y values
init_matrix = [];          % matrix of initial points

for ci=1:length(initXLim)
    for cj=1:length(initYLim)
        init = [initXLim(ci) initYLim(cj)];
        init_matrix = [init_matrix; init]; 
    end
end

%optimization
xmin  = zeros(size(init_matrix));

for ck=1:length(init_matrix)
    xmin(ck,:) = fminunc(F, init_matrix(ck,:), options);
end

% find the convergence points of initial values
xinit_fix1_dist = sqrt((xmin(:,1)-xmin1(1)).^2 + (xmin(:,2)-xmin1(2)).^2); % distance to fix1 at each fix
xinit_fix2_dist = sqrt((xmin(:,1)-xmin2(1)).^2 + (xmin(:,2)-xmin2(2)).^2); % distance to fix2 at each fix
xinit_fix1      = xinit_fix1_dist <= xinit_fix2_dist;                      % converges to fix1 if closer to it
xinit_fix2      = xinit_fix1_dist >  xinit_fix2_dist;                      % converges to fix2 if closer to it
xM              = init_matrix(xinit_fix1,:);                               % initial values that converge to fix1
xC              = init_matrix(xinit_fix2,:);                               % initial values that converge to fix2

% plots
figure(3)
viscircles(L1,R1,'Color','g');
viscircles(L2,R2,'Color','r');
viscircles(L3,R3,'Color','b');
hold on
scatter(xmin1(1),xmin1(2),40,'m','fill')
scatter(xmin2(1),xmin2(2),40,'c','fill')
% scatter(xmin3(1),xmin3(2),40,'k','fill')
title('Problem 3: Part c','Inconsistent Case: Fixes From a Grid of Initial Points')
xlabel('x')
ylabel('y')
legend('Fix 1','Fix 2','Location','best')
grid on

figure(4)
scatter(xM(:,1),xM(:,2),'m','fill')
hold on
scatter(xC(:,1),xC(:,2),'c','fill')
title('Problem 3: Part c','Inconsistent Case: Initial Points')
xlabel('x_{init}')
ylabel('y_{init}')
legend('Fix 1','Fix 2','Location','best')
grid on


