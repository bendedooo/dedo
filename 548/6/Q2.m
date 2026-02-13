clear all; clc; close all;
rho = 0.01215;
x = 0.75; y = -0.05; z = 0;
r1 = sqrt((x+rho)^2+y^2+z^2);
r2 = sqrt((x-1+rho)^2+y^2+z^2);
U = -0.5*((x^2+y^2)+ 2*(1-rho)/r1 + 2*rho/r2 + rho*(1-rho))
U1 = -1.60017;
v = sqrt(2*(U1-U))