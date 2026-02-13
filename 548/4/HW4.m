%% Problem 1

% Part c)

clc; clear all; close all;
mu = 398600.4405;
V_inf = 2.98;
R_0 = 384000;
R_p = 6578;

deltaVa = sqrt(V_inf^2+2*mu/R_0)-sqrt(mu/R_0)
deltaV1 = sqrt(mu/R_0) - sqrt(2*mu*(1/R_0 - 1/(R_0+R_p)));
deltaV2 = sqrt(V_inf^2+2*mu/R_p) - sqrt(2*mu*(1/R_p - 1/(R_0+R_p)));
deltaVb = deltaV1+deltaV2

% part d
R_0 = 10000;
deltaVa = sqrt(V_inf^2+2*mu/R_0)-sqrt(mu/R_0)
deltaV1 = sqrt(mu/R_0) - sqrt(2*mu*(1/R_0 - 1/(R_0+R_p)));
deltaV2 = sqrt(V_inf^2+2*mu/R_p) - sqrt(2*mu*(1/R_p - 1/(R_0+R_p)));
deltaVb = deltaV1+deltaV2


%% PROBLEM 2
Re      = 6378;      % Earth Radius (km)
mu      = 398600.4405;  % Gravitational constant for Earth (km^3/s^2) 
m      = 100; %mass [kg]
% Initial position vector in ECI frame (units: km)
R0 = [6628; 0; 0]; % [-8.058685464116713e+002;-4.135050989428283e+003;5.071447347723579e+003];   
% Initial velocity vector in ECI frame (units: km/sec)
V0 = [0; 5.4836; 5.4836]; 
% Initial state [position vector on top of velocity vector] 
X0 = [R0; V0]; 

% For spherical coord
r_S = norm(R0);
nu_S = atan2(R0(2),R0(1));
lambda_S = atan2( R0(3),(sqrt(R0(1)^2+(R0(2)^2))) );
RS0   = [r_S; nu_S; lambda_S];
C_G2S = [ cos(nu_S)*cos(lambda_S)    sin(nu_S)*cos(lambda_S)   sin(lambda_S);
         -sin(nu_S)                  cos(nu_S)                 0;
         -cos(nu_S)*sin(lambda_S)   -sin(nu_S)*sin(lambda_S)   cos(lambda_S)];
VS0   = C_G2S*V0;
RS0dot(1) = VS0(1);
RS0dot(2) = VS0(2)/cos(RS0(3))/RS0(1);
RS0dot(3) = VS0(3)/RS0(1);
% Initial states in spherical coord
XS0   = [RS0; RS0dot'];

t0     = 0;          % initial time, sec
tstep  = 20;         % time step to sample the output, sec
tfinal = 108000;     % final time, sec, approximately 5 orbits
 
% Integrate ODEs using Runge-Kutta 45 method
tspan    = t0:tstep:tfinal;
options = odeset('AbsTol',1e-6,'RelTol',1e-8); 
[Tp, Xp] = ode45(@satdyn,tspan,X0,options);
timeVec = Tp;XTraj = Xp;

% post process
r_S_ECI = []; nu_S_ECI = []; lambda_S_ECI = [];
for k=1:1:length(timeVec)

    Xk = XTraj(k,:);
    r_S_ECI(k) = norm(Xk(1:3));
    nu_S_ECI(k) = atan2(Xk(2),Xk(1));
    lambda_S_ECI(k) = atan2( Xk(3),(sqrt(Xk(1)^2+(Xk(2)^2))) );

end

% Integrate ODEs using Runge-Kutta 45 method SPHERICAL COORD
tspan    = t0:tstep:tfinal;
options = odeset('AbsTol',1e-6,'RelTol',1e-8); 
[TSp, XSp] = ode45(@satdynS,tspan,XS0,options);
timeVecS = TSp;XTrajS = XSp;

% post process
x_S = []; y_S = []; z_S = [];
for k=1:1:length(timeVecS)
    Xk = XTrajS(k,:);
    rk = Xk(1);
    nuk = Xk(2);
    lambdak = Xk(3);
    x_S(k) = rk*cos(nuk)*cos(lambdak);
    y_S(k) = rk*sin(nuk)*cos(lambdak);
    z_S(k) = rk*sin(lambdak);
end

% plots
figure(1);
plot3(Xp(:,1), Xp(:,2), Xp(:,3),'g-','linewidth',3); 
xlabel('x (km)','fontsize',12);ylabel('y (km)','fontsize',12); 
zlabel('z (km)','fontsize',12); 
set(gca,'fontsize',12)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
title('3D Trajectory Cartesian Coordinates')
hold off

% plot altitude above the Earth, km
figure(2);
plot(Tp/60, sqrt(Xp(:,1).^2 + Xp(:,2).^2 + Xp(:,3).^2) - Re, 'linewidth',3);
ylabel('Altitude r - R_e (km)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Altitude r - R_e (km) Cartesian Coordinates')

figure(3);
plot(TSp/60, wrapToPi(nu_S_ECI), 'linewidth',3);
ylabel('Right Ascention (rad)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Right Ascention (rad) ECI')

figure(4);
plot(TSp/60, lambda_S_ECI, 'linewidth',3);
ylabel('Declination (rad)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Declination (rad) ECI')

figure(5);
plot3(x_S, y_S, z_S,'g-','linewidth',3); 
xlabel('x (km)','fontsize',12);ylabel('y (km)','fontsize',12); 
zlabel('z (km)','fontsize',12); 
set(gca,'fontsize',12)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
title('3D Trajectory Spherical Coordinates')
hold off

figure(6);
plot(TSp/60, XSp(:,1)- Re, 'linewidth',3);
ylabel('Altitude r - R_e (km)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Altitude r - R_e (km) Spherical Coordinates')

figure(7);
plot(TSp/60, wrapToPi(XSp(:,2)), 'linewidth',3);
ylabel('Right Ascention (rad)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Right Ascention (rad) Spherical Coordinates')

figure(8);
plot(TSp/60, XSp(:,3), 'linewidth',3);
ylabel('Declination (rad)','fontsize',12)
xlabel('t (min)','fontsize',12)
set(gca,'fontsize',12)
title('Declination (rad) Spherical Coordinates')

function Xdot = satdyn(t,X) 
mu = 398600.4405;        % gravitational constant for Earth (km^3/s^2) 
r  = X(1:3);          % position vector (km)
v  = X(4:6);          % velocity vector (km/sec) 
h  = cross(r,v);
x = r(1);
y = r(2);
z = r(3);
F_p = 0.005*v/norm(v);
rn  = norm(r); % distance (km)
m = 100; %mass [kg]
Xdot(1:3)  = v;
Xdot(4:6)  = -mu/rn^3*r + F_p/m; 
Xdot       = Xdot(:); % convert to vector column format

end

function Xdot = satdynS(t,X) 
mu = 398600;%.4405;        % gravitational constant for Earth (km^3/s^2) 
m = 100; %mass [kg]

V_r = X(4);
V_nu = X(1)*X(5)*cos(X(3));
V_lambda = X(1)*X(6);
V_S = [V_r; V_nu; V_lambda];
e_v = V_S/norm(V_S);
F_p = 0.005*e_v/m;

Xdot(1:3)  = X(4:6); 
% Xdot(4:6)  = -mu/X(1)^2.*[1; 0; 0] + F_p/m; 
% Xdot       = Xdot(:); % convert to vector column format
r = X(1);
nu = X(2);
% nu = wrapToPi(nu);
lambda = X(3);
r_dot = X(4);
nu_dot = X(5);
lambda_dot = X(6);

F_r = F_p(1);
F_nu = F_p(2);
F_lambda = F_p(3);

Xdot(4) = -mu/r^2 + r*nu_dot^2*cos(lambda)^2 + r*lambda_dot^2 + F_r;
Xdot(5) = (r*nu_dot*sin(lambda)*lambda_dot +r*lambda_dot*nu_dot*sin(lambda) - 2*r_dot*nu_dot*cos(lambda) + F_nu)/(r*cos(lambda));
Xdot(6) = (-2*r_dot*lambda_dot - r*nu_dot^2*cos(lambda)*sin(lambda) + F_lambda)/r;

Xdot       = Xdot(:); % convert to vector column format
end
