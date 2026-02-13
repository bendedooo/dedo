close all; clear all;
global rho
rho = 1.215059e-2; % Earth-Moon CR3BP
global part
part = 'c'; % 'a', 'b' or 'c'

% Note state order: [x, y, z, xdot, ydot ,zdot]

tfinal = 2.7462016488; % L1 Halo Orbit
if part == 'c'
    tfinal = tfinal*4;
end
X0     = [0.8368126154 0 -0.1474695518 0  0.2560040701 0]';

% Parameters for integrating 3BP ODEs using Runge-Kutta 45 method
tstep  = 0.001;   % time step to sample the output, scaled units
tspan = [0: tstep: tfinal];
options = odeset('AbsTol',1e-11,'RelTol',1e-13);

%% Simulate actual and reference trajectories and plot required quantities

% Initialization X_init = [X_r(0); X_a(0)];
X_r0 = X0;
X_a0 = X_r0 + [-0.02 -0.02 0 0 0 0]';

% Use ode45 and 'ode3BPclp' function to simulate reference trajectory
% X_r(t) and X_a(t)
[TOut, XOut] = ode45(@ode3BPclp,tspan,[X_r0; X_a0],options);

% Post-processing and plotting results
D = 384748; %from module 7 page 8
XOut(:,1:3) = XOut(:,1:3) * D;
XOut(:,7:9) = XOut(:,7:9) * D;

% n = 2.661699e-06; %from module 7 page 8
% XOut(:,4:6) = XOut(:,4:6) * n*D;
% XOut(:,10:12) = XOut(:,10:12) * n*D;

figure(1);
plot3(XOut(:,1), XOut(:,2), XOut(:,3),'g-','linewidth',3);
hold on
plot3(XOut(:,7), XOut(:,8), XOut(:,9),'b-.','linewidth',3);
plot3(XOut(1,1), XOut(1,2), XOut(1,3),'g*','linewidth',3);
plot3(XOut(1,7), XOut(1,8), XOut(1,9),'b*','linewidth',3);
xlabel('x (km)','fontsize',12);
ylabel('y (km)','fontsize',12);
zlabel('z (km)','fontsize',12);
% set(gca,'fontsize',12)
% axis equal
title('3D Trajectories')
legend('Reference Trajectory','Actual Trajectory', Location='best');%,'Reference Trajectory Initial Point','Actual Trajectory Initial Point')
%% Implement 'ode3BPclp'

% ODE file for 3BP. NOTE: Time and Distance are SCALED
% 1 unit of time = 1/n sec
% 1 unit of distance = D km, D=distance between primaries
% input:    scaled time t, augmented state vector X = [X_r(t); X_a(t)];
%           X_r(t) = [x_r(t), y_r(t), z_r(t), xrdot(t), yrdot(t), zrdot(t)]
%           X_a(t) = [x_a(t), y_a(t), z_a(t), xadot(t), yadot(t), zadot(t)]
% output:   Xdot = [Xrdot(t); Xadot(t)];
%           ux, uy, uz: dimensionless force Fx, Fy, Fz divided by m*n^2*D
function [Xdot, ux, uy, uz] = ode3BPclp(t, X)

global rho
global part

Xr = X(1: 6); Xa = X(7: 12);

% compute Xrdot(t) based on the unforced EoMs
%add your code and modify the next line%
Xrdot = CR3BP(Xr,[0 0 0]);

% compute Xadot(t) by the following steps
% generate ux, uy, and uz based on feedback law
if part == 'a' || part == 'c'
    K = [9.9504 , -2.8108 , 1.3799 , 4.3845 , 0.5600 , 0.3518 ;
        3.8239 , 2.0919 , 0.3655 , 0.5600 , 2.8454 , 0.1188 ;
        1.5558 , -0.1666 , 3.8614 , 0.3518 , 0.1188 , 3.0959];
    u = -K*(Xa-Xr);

elseif part == 'b'
    u = [0 0 0]';
end
ux = u(1);
uy = u(2);
uz = u(3);

% compute time rate of change of the spacecraft state, modify the next line%
Xadot = CR3BP(Xa,u);

% compute Xdot
Xdot = [Xrdot(:); Xadot(:)];

end

%% Implement ode function 'CR3BP'
% Dimensionless ode for CR3BP problem in rotating frame and in scaled units

function Xdot = CR3BP(X, u)

global rho

x = X(1); y = X(2); z = X(3); xdot = X(4); ydot = X(5); zdot = X(6);
ux = u(1); uy = u(2); uz = u(3);
r1 = sqrt((x+rho)^2+y^2+z^2);
r2 = sqrt((x-1+rho)^2+y^2+z^2);

xddot = 2*ydot + x -(1-rho)*(x+rho)/r1^3 -(rho)*(x-1+rho)/r2^3 + ux;
yddot = -2*xdot + y -(1-rho)*(y)/r1^3 -(rho)*(y)/r2^3 + uy;
zddot = -(1-rho)*(z)/r1^3 -(rho)*(z)/r2^3 + uz;

Xdot = [xdot,ydot,zdot,xddot,yddot,zddot]';
end
