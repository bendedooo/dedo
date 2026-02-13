clear all; clc; close all;
rho = 0.953875e-03;

fun = @LagrangePoints;
x0 = 0;
xe1 = fsolve(fun,x0);
L1 = [xe1, 0, 0];

x0 = 1;
xe2 = fsolve(fun,x0);
L2 = [xe2, 0, 0];

x0 = -1.1984;
xe3 = fsolve(fun,x0);
L3 = [xe3, 0, 0];

L4 = [0.5-rho, sqrt(3)/2, 0];
L5 = [0.5-rho, -sqrt(3)/2, 0];

L = [L1;L2;L3;L4;L5]

for i=1:5
    Lpt = L(i,:);
    CJ = JacobiConst(Lpt)
end

function F = LagrangePoints(x)
rho = 0.953875e-03;
F(1) = x - (1-rho)*(x+rho)/(abs(x+rho))^3 - rho*(x-1+rho)/(abs(x-1+rho))^3;
end

function F = JacobiConst(L)
rho = 0.953875e-03;
x=L(1);
y=L(2);
z=L(3);
r1 = sqrt((x+rho)^2+y^2+z^2);
r2 = sqrt((x-1+rho)^2+y^2+z^2);
F(1) = (x^2+y^2)+ 2*(1-rho)/r1 + 2*rho/r2 + rho*(1-rho);
end