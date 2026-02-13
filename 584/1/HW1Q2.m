% PROBLEM 2
clc
clear all
close all
format long

% beacons and measurements
L1 = [2 2];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
R1 = 4.8;    % range measurement from L1

L2 = [6 6];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
R2 = 8.4;      % range measurement from L2

% initial position guesses
init1 = [14 2];
init2 = [-2 16];

global h;    % array to store position fixes at each iteration
h = [];

F       = @(x) HW1Q2costFun(x,L1,L2,R1,R2);                   % define the cost function
%% PART A

% optimization
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12); % set tolerance

[xmin1, fval1] = fminunc(F, init1, options);                  
[xmin2, fval2] = fminunc(F, init2, options);

% plot
figure(1)
viscircles(L1,R1,'Color','g');
viscircles(L2,R2,'Color','r');
hold on
scatter(xmin1(1),xmin1(2),40,'m','fill')
scatter(xmin2(1),xmin2(2),40,'c','fill')
scatter(init1(1),init1(2),40,'m*')
scatter(init2(1),init2(2),40,'c*')
title('Problem 2: Part a','2D Position Fix with Two Range Measurements')
xlabel('x')
ylabel('y')
legend('Fix 1','Fix 2','Initial Guess 1','Initial Guess 2','Location','best')
grid on

% [xout,yout] = circcirc(L1(1),L1(2),R1,L2(1),L2(2),R2); % intersection point(s) of two circles

%% PART B

% optimization
options      = optimoptions(@fminunc,'OptimalityTolerance',1e-12,'OutputFcn',@HW1Q2outfun);  %  set the tolerance and outfun to store iteration outputs
[xmin, fval] = fminunc(F, init2, options);

% plot iteration vs distance
iter = 1:length(h(:,1));  % number of iterations

figure(2)
semilogy(iter, sqrt((h(:,1)-xmin(1)).^2 + (h(:,2)-xmin(2)).^2),'b','linewidth',2); 
xlim([min(iter) max(iter)]);
xticks(1:2:max(iter)); 
title('Problem 2: Part b','Iteration vs Distance Error')
xlabel('Iteration')
ylabel('Distance Error')
grid on

%% PART C 

% crate a grid of initial points
initXLim    = [-8:2:12];   % range for initial point x values
initYLim    = [-8:2:12];   % range for initial point y values
init_matrix = [];          % matrix of initial points

for ci=1:length(initXLim)
    for cj=1:length(initYLim)
        init = [initXLim(ci) initYLim(cj)];
        init_matrix = [init_matrix; init]; 
    end
end

% optimization
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12);    % set the tolerance
xmin    = zeros(size(init_matrix));

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
hold on
scatter(xmin1(1),xmin1(2),40,'m','fill')
scatter(xmin2(1),xmin2(2),40,'c','fill')
title('Problem 2: Part c','Fixes From a Grid of Initial Points')
xlabel('x')
ylabel('y')
legend('Fix 1','Fix 2','Location','best')
grid on

figure(4)
scatter(xM(:,1),xM(:,2),'m','fill')
hold on
scatter(xC(:,1),xC(:,2),'c','fill')
title('Problem 2: Part c','Initial Points')
xlabel('x_{init}')
ylabel('y_{init}')
legend('Fix 1','Fix 2','Location','best')
grid on