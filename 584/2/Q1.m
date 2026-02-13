% Q1

clc
clear all
close all
format long

rad = pi/180;

%% PART C
% beacons and measurements
L1 = [0 0];  % first known location coordinates [x y]
x1 = L1(1);
y1 = L1(2);
theta1 = -165*rad;    % measurement from L1

L2 = [4 2];  % second known location coordinates [x y]
x2 = L2(1);
y2 = L2(2);
theta2 = 150*rad;      % measurement from L2

T1 = tan(pi/2-theta1);
T2 = tan(pi/2-theta2);
A=[T1 -1;T2 -1];
B=[T1*x1-y1;T2*x2-y2];
C=inv(A)*B;
x0=C(1);
y0=C(2);


% plot
figure(1)
scatter(x1,y1,40,'g','fill')
hold on
scatter(x2,y2,40,'g','fill')
scatter(x0,y0,40,'r','fill')
title('Problem 1: Part c','2D Position Fix with Two Beacon Measurements')
xlabel('x [m]')
ylabel('y [m]')
legend('Lighthouse 1','Lighthouse 2','Position Fix','Location','best')
grid on

%% PART D

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

T1 = tan(pi/2-theta1);
T2 = tan(pi/2-theta2);
T3 = tan(pi/2-theta3);
T  = [T1 T2 T3];
x0y0 = [];
x = [x1 x2 x3];
y = [y1 y2 y3];
for i=1:length(T)
    for j=1:length(T)
        if i<j || i==~j
A=[T(i) -1; T(j) -1];
B=[T(i)*x(i)-y(i);T(j)*x(j)-y(j)];
C=inv(A)*B;
x0y0 = [x0y0; C(1) C(2)];
        end
    end
end

% plot
figure(2)
scatter(x,y,40,'g','fill')
hold on
scatter(x0y0(:,1),x0y0(:,2),40,'b','fill')
scatter(sum(x0y0(:,1))/length(x0y0(:,1)),sum(x0y0(:,2))/length(x0y0(:,2)),40,'r','fill')
title('Problem 1: Part d','3 Lighthouses')
xlabel('x [m]')
ylabel('y [m]')
legend('Lighthouses','Position Fixes' ,'Lighthouse Center','Location','best')
grid on