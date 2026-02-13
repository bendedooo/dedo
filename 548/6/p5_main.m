close all; clear all;

global rho
rho = 1.215059e-2; % Earth-Moon CR3BP

% Note state order: [x, y, z, xdot, ydot ,zdot]

tfinal = 2.7462016488; % L1 Halo Orbit
X0     = [0.8368126154 0 -0.1474695518 0  0.2560040701 0]';

% Parameters for integrating 3BP ODEs using Runge-Kutta 45 method
tstep  = 0.001;   % time step to sample the output, scaled units
tspan = [0: tstep: tfinal];
options = odeset('AbsTol',1e-11,'RelTol',1e-13); 

%% Simulate actual and reference trajectories and plot required quantities

% Initialization X_init = [X_r(0); X_a(0)];

% Use ode45 and 'ode3BPclp' function to simulate reference trajectory
% X_r(t) and X_a(t)

% Post-processing and plotting results


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

Xr = X(1: 6); Xa = X(7: 12);

% compute Xrdot(t) based on the unforced EoMs 
    %add your code and modify the next line%
    Xrdot = CR3BP([],[]);

% compute Xadot(t) by the following steps

    % generate ux, uy, and uz based on feedback law

    % compute time rate of change of the spacecraft state, modify the next line%
    Xadot = CR3BP([],[]);

% compute Xdot
    Xdot = [Xrdot(:); Xadot(:)];

end

%% Implement ode function 'CR3BP'
% Dimensionless ode for CR3BP problem in rotating frame and in scaled units

function Xdot = CR3BP(X, u)

global rho

x = X(1); y = X(2); z = X(3); xdot = X(4); ydot = X(5); zdot = X(6);
ux = u(1); uy = u(2); uz = u(3);

% ADD YOUR CODE HERE %

Xdot = [xdot,ydot,zdot,xddot,yddot,zddot]';
end
