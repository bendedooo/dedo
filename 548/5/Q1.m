%% PROBLEM 1
clc
clear all
close all

Re      = 6378;      % Earth Radius (km)
mu      = 398600.4405;  % Gravitational constant for Earth (km^3/s^2) 
m = 100; %mass [kg]
% Initial position vector in ECI frame (units: km)
R0 = [6628; 0; 0]; % [-8.058685464116713e+002;-4.135050989428283e+003;5.071447347723579e+003];  
% Initial velocity vector in ECI frame (units: km/sec)
V0 = [0; 5.4836; 5.4836]; 
 
% Initial state [position vector on top of velocity vector] 
X0 = [R0; V0]; 

% Cylindrical coordinates
i_rho = [1 0 0];
i_phi = [0 1 0];
i_z   = [0 0 1];

rho = dot(R0,i_rho);
phi = 0;
z   = dot(R0,i_z);

R0_C =[rho; phi; z];

rho_dot = dot(V0,i_rho);
phi_dot = dot(V0,i_phi) / rho;
z_dot   = dot(V0,i_z);

V0_C    = [rho_dot; phi_dot; z_dot];
% state vector in cylindrical coords
X0_C = [R0_C;V0_C];

t0     = 0;          % initial time, sec
tstep  = 20;         % time step to sample the output, sec
tfinal = 108000;     % final time, sec, approximately 5 orbits
 
% Integrate ODEs using Runge-Kutta 45 method
tspan    = t0:tstep:tfinal;
options = odeset('AbsTol',1e-6,'RelTol',1e-8); 
[Tp, Xp] = ode45(@satdyn,tspan,X0_C,options);
timeVec = Tp;XTraj = Xp;

figure(1)
plot(Tp/60, Xp(:,1), 'linewidth',3);
ylabel('\rho','fontsize',12)

figure(2)
plot(Tp/60, wrapTo2Pi(Xp(:,2)), 'linewidth',3);
ylabel('\phi','fontsize',12)

figure(3)
plot(Tp/60, Xp(:,3), 'linewidth',3);
ylabel('z','fontsize',12)

x_ECI=[];y_ECI=[];z_ECI=[];
for k=1:1:length(timeVec)

    Xk = XTraj(k,:);
    x_ECI(k) = Xk(1)*cos(Xk(2));
    y_ECI(k) = Xk(1)*sin(Xk(2));
    z_ECI(k) = Xk(3);

end

figure(4);
plot3(x_ECI, y_ECI, z_ECI,'g-','linewidth',3); 
xlabel('x (km)','fontsize',12);ylabel('y (km)','fontsize',12); 
zlabel('z (km)','fontsize',12); 
set(gca,'fontsize',12)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
title('3D Trajectory')
hold off

% ylabel('Declination (rad)','fontsize',12)
% xlabel('t (min)','fontsize',12)
% set(gca,'fontsize',12)
% title('Declination (rad) Spherical Coordinates')



function Xdot = satdyn(t,X) 
mu = 398600;%.4405;        % gravitational constant for Earth (km^3/s^2) 
m = 100; %mass [kg]

rho = X(1);
phi = X(2);
z   = X(3);
rhodot = X(4);
phidot = X(5);
zdot   = X(6);

i_rho =  cos(phi)*[1 0 0] + sin(phi)*[0 1 0];
i_phi = -sin(phi)*[1 0 0] + cos(phi)*[0 1 0];
i_z   = [0 0 1];

V_S = rhodot*i_rho+ rho*phidot*i_phi +zdot*i_z;
% V_S = X(4:6).*[1;rho;1];
e_v = V_S/norm(V_S);


F_p = 0.005*e_v/m;
 
F_rho = dot(F_p,i_rho);
F_phi = dot(F_p,i_phi);
F_z   = dot(F_p,i_z);


% V_S = [rhodot; rho*phidot; zdot];
% e_v = V_S/norm(V_S);
% F_p = 0.005*e_v/m;
%  
% F_rho = F_p(1);
% F_phi = F_p(2);
% F_z   = F_p(3);


Xdot(1:3)  = X(4:6);
Xdot(4) = F_rho+rho*phidot^2-mu*rho/(rho^2+z^2)^(1.5);
Xdot(5) = (F_phi*rho-2*rho*rhodot*phidot)/(rho^2);
Xdot(6) = F_z-mu*z/(rho^2+z^2)^(1.5);

Xdot       = Xdot(:); % convert to vector column format

end





