clear all; clc; close all
% Question 5 Part c
d2r = pi/180;
T_max = 11000; %[lbs]
W = 75000; %[lbs]
rho = 0.00238; %[slug/ft^3]
C_D_0 = 0.015;
S = 950; %[ft^2]
C_L_alpha = 0.12;
K = 0.05;

x0 = [0;1000;0;0;0]; % gamma v alpha p1 p2
% 
% F = @(alpha,V,gamma) [-gamma + p1*dg1_alpha + p2*dg2_alpha;...
%                       -gamma + p1*dg1_V + p2*dg2_V;...
%                       -gamma + p1*dg1_gamma + p2*dg2_gamma;...
%                       g1;...
%                       g2];

F = @(x) [-1 + x(4)*(-(1250*pi*cos((pi*x(1))/180))/3)                                 + x(5)*(1250*pi*sin((pi*x(1))/180))/3;...
           0 + x(4)*(-(2261*x(2)*((9*x(3)^2)/12500 + 3/200))/1000)                    + x(5)*(6783*x(2)*x(3))/25000;...
           0 + x(4)*(- (550*pi*sin((pi*x(3))/180))/9 - (20349*x(2)^2*x(3))/12500000)  + x(5)*((6783*x(2)^2)/50000 + (550*pi*cos((x(3)*pi)/180))/9);...
          T_max*cos(x(3)*d2r) - 0.5*rho*x(2)^2*S*(C_D_0 + K*(C_L_alpha*x(3))^2) - W*sin(x(1)*d2r);...
          0.5*rho*x(2)^2*S*C_L_alpha*x(3) + T_max*sin(x(3)*d2r)- W*cos(x(1)*d2r)];

options=optimset('MaxFunEvals',100000,'MaxIter', 100000,'TolFun',1e-10,'TolX',1e-6, 'Display','final-detailed');
% options = optimoptions('fsolve','Display','iter');

[X,fval] = fsolve(F,x0,options)


