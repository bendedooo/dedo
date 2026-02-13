% PROBLEM 4
clc
clear all
close all
format long

% beacons and measurements same as question 3
load Q3
L1     = L1;    % first known location coordinates [x y]
x1     = L1(1);
y1     = L1(2);
R1true = R1;       % range measurement from L1

L2     = L2;    % second known location coordinates [x y]
x2     = L2(1);
y2     = L2(2);
R2true = R2;        % range measurement from L2

L3     = L3;  % third known location coordinates [x y]
x3     = L3(1);
y3     = L3(2);
R3true = R3;       % range measurement from L3

% initial position guess
init1 = [14 2]; 

%% PART A 

a = 3; % noise coefficient
xmin    = zeros(100,2);

for i=1:100
% generate noised measurements
R1 = R1true + a*(rand(1) - 0.5);
R2 = R2true + a*(rand(1) - 0.5);
R3 = R3true + a*(rand(1) - 0.5);

%optimization
F = @(x) HW1Q3costFun(x,L1,L2,L3,R1,R2,R3);                   % new cost function
options = optimoptions(@fminunc,'OptimalityTolerance',1e-12); % set the tolerance

[xmin(i,:), fval1] = fminunc(F, init1, options);
end

figure(1)
viscircles(L1,R1true,'Color','g');
viscircles(L2,R2true,'Color','r');
viscircles(L3,R3true,'Color','b');
hold on
scatter(xmin(:,1),xmin(:,2),40,'r','fill')
scatter(xmin1(:,1),xmin1(:,2),40,'m','fill')
title('Problem 4: Part a','2D Position Fix with Three Noisy Range Measurements')
xlabel('x')
ylabel('y')
legend('Fix_i , i=1,...,100','True Fix','Location','best')
grid on

%% PART B 

xmin      = zeros(100,2);
dist2true = [];        % distance to true fix

a = [0:0.15:20];   % noise coefficients
randmatrix = [rand(1,100)' rand(1,100)' rand(1,100)']; % noise random number matrix


for ib=1:length(a)
        for j=1:100;
        R1 = R1true + a(ib)*(randmatrix(j,1) - 0.5);
        R2 = R2true + a(ib)*(randmatrix(j,2) - 0.5);
        R3 = R3true + a(ib)*(randmatrix(j,3) - 0.5);
        % optimization
        F = @(x) HW1Q3costFun(x,L1,L2,L3,R1,R2,R3);
        options = optimoptions(@fminunc,'OptimalityTolerance',1e-12);                % set the tolerance
        [xmin(j,:), fval1] = fminunc(F, init1, options);
        fixnoised_dist(j) = sqrt((xmin(j,1)-xmin1(1)).^2 + (xmin(j,2)-xmin1(2)).^2); % distance to true fix
        end
    dist2true = [dist2true fixnoised_dist']; % distance to true fix at noised measurement with a
end

[maxDist,indMaxDist] = max(dist2true);       % [maximum distance to true fix for each a, index for maxDist]

figure(2)
% semilogy(a,maxDist,'r*')
semilogy(a,maxDist,'LineWidth',2)
title('Problem 4: Part b','Max Distance to the True Fix for 2D Position Fix with Different Noisy Range Measurements')
xlabel('a value')
ylabel('Distance_{max} to True Fix')
legend('a_i','Location','best')
grid on

