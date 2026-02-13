clear all; clc; close all;
%% part b
fun = @(x) -5*x(1) - x(2);
x0 = [0,0];
A = [-1  0;...
      3  1;...
      1 -2];
b = [0; 11; 2];
x = fmincon(fun,x0,A,b)

%% part c
f = [-5 -1];
x = linprog(f,A,b)

%% part e

tol = 1e-01;
actConst = abs(A*x-b) < tol
isFull = 2 == rank([3 1;1 -2])
