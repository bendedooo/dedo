clear all; clc; close all;

x0N = [-10; -10];
x0P = [ 10;  10];

A = [];
b = [];
Aeq = [];
beq = [];

lb = [];
ub = [];

options=optimset('MaxFunEvals',100000,'MaxIter', 100000,'TolFun',1e-10,'TolX',1e-6, 'Display','final-detailed');

x0N
xN = fmincon(@fun,x0N,A,b,Aeq,beq,lb,ub,@mycon,options)
x0P
xP = fmincon(@fun,x0P,A,b,Aeq,beq,lb,ub,@mycon,options)


function out = fun(x)
out = x(1)^2 + x(2)^2;
end

function [c,ceq] = mycon(x)
c = 1 - x(1)*x(2);     % Compute nonlinear inequalities at x.
ceq = [];              % Compute nonlinear equalities at x.
end
