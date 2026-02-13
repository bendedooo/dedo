clear all; clc; close all;
% Question 5 part a
% initial values
gamma0 = 0;
V0     = 1000;
alpha0 = 0;

x0 = [gamma0 V0 alpha0];

A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [];
ub  = [];
% lb = [-180 -inf -180];
% ub = [180 inf 180];

options=optimset('MaxFunEvals',100000,'MaxIter', 100000,'TolFun',1e-10,'TolX',1e-6, 'Display','final-detailed');
x = fmincon(@fun,x0,A,b,Aeq,beq,lb,ub,@mycon,options)

function out = fun(x)
% gamma, V, alpha = [x(1) x(2) x(3)];
gamma = x(1);
V     = x(2);
alpha = x(3);

% output
out   = -gamma;
% out2 = -x(1);
end

function [c,ceq] = mycon(x)

d2r       = pi/180;
T_max     = 11000; %[lbs]
W         = 75000; %[lbs]
rho       = 0.00238; %[slug/ft^3]
C_D_0     = 0.015;
S         = 950; %[ft^2]
C_L_alpha = 0.12;
K         = 0.05;

gamma = x(1);
V     = x(2);
alpha = x(3);

alpha_r = alpha*d2r;
gamma_r = gamma*d2r;

C_L = C_L_alpha*alpha;
C_D = C_D_0 + K*C_L^2;

D = 0.5*rho*V^2*S*C_D;
L = 0.5*rho*V^2*S*C_L;

% outputs
c = [];     % Compute nonlinear inequalities at x.
ceq(1) = T_max*cos(alpha_r) - D - W*sin(gamma_r);              % Compute nonlinear equalities at x.
ceq(2) = L + T_max*sin(alpha_r)- W*cos(gamma_r);

end


