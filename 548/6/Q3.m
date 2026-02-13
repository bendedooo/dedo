% HW6Q3
clear all; clc; close all;
rho = 0.03;

% Initial states 
X0   = [0.46;sqrt(3)/2;0;0]; % states x,y,xdot,ydot

t0     = 0;          % initial time, sec
tstep  = 0.2;         % time step to sample the output, sec
tfinal = 60*pi;     % final time, sec, approximately 5 orbits

tspan    = t0:tstep:tfinal;
options = odeset('AbsTol',1e-6,'RelTol',1e-8); 
[Tp, Xp] = ode45(@CR3BP,tspan,X0,options);
timeVec = Tp;XTraj = Xp;

% plots
figure(1);
plot(Xp(:,1), Xp(:,2), 'g-','linewidth',3); 
xlabel('x ','fontsize',12);
ylabel('y ','fontsize',12); 
set(gca,'fontsize',12)
hold on;
plot(rho, 0, 'b*','linewidth',3); 
plot(1-rho,0, 'r*','linewidth',3); 
axis equal
legend('Trajectory','Fixed Body 1','Fixed Body 2')
title('2D Trajectory')
hold off

function Xdot = CR3BP(t,X) 

rho = 0.03;
x = X(1);
y = X(2);
z = 0;

r1 = sqrt((x+rho)^2+y^2+z^2);
r2 = sqrt((x-1+rho)^2+y^2+z^2);

Xdot(1:2)  = X(3:4); 
Xdot(3)    = 2*Xdot(2)  + x - (1-rho)*(x+rho)/r1^3 - rho*(x-1+rho)/r2^3;
Xdot(4)    = -2*Xdot(1) + y*(1 - (1-rho)/r1^3 - rho/r2^3);

Xdot       = Xdot(:); % convert to vector column format
end