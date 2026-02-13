%% problem 1
clc; clear all; close all;
mu = 398600.4405;
thetadot = 7.2912e-05;
phi = pi/4;
r = 42189.7;

e = norm((thetadot*r*cos(phi))^2*r/mu -1)

p = (thetadot*cos(phi)*r^2)^2/mu

rp = p/(1+e)
ra = p/(1-e)
a = (ra+rp)/2

%% problem 3
clc; clear all; close all;
Re  = 6378;
mu  = 398600.4405;
ha  = 400:120000;
hp  = 400;
Rp  = hp + Re;
Ra  = ha + Re;
ecc = (Ra - Rp)./(Ra + Rp);
a   = (Ra + Rp)./2;
n   = sqrt(mu ./ (a.^3));
rangeDeg = 20;
theta1   = (180-rangeDeg)*pi/180;
theta2   = (180+rangeDeg)*pi/180;
eccAnom1 = 2*atan(tan(theta1/2)*sqrt((1-ecc)./(1+ecc)));
eccAnom2 = 2*pi + 2*atan(tan(theta2/2)*sqrt((1-ecc)./(1+ecc)));
meanAnom1 = eccAnom1 - ecc.*sin(eccAnom1);
meanAnom2 = eccAnom2 - ecc.*sin(eccAnom2);
TIM = (meanAnom2-meanAnom1)./n; %time in orbit
T = 2*pi./n; %period
timeFrac = TIM./T;

figure(1);
plot(ha, timeFrac, 'linewidth',3);
ylabel('Time of Flight Fraction','fontsize',14)
xlabel('h_a (km)','fontsize',14)
title('Time of Flight Fraction vs h_a')
set(gca,'fontsize',14)

figure(2);
plot(ecc, timeFrac, 'linewidth',3);
ylabel('Time of Flight Fraction','fontsize',14)
xlabel('Eccentricity','fontsize',14)
title('Time of Flight Fraction vs e')
set(gca,'fontsize',14)

%% problem 4
clc; clear all; close all;
% part a
syms r x y z R mu J2 J3
r = sqrt(x^2+y^2+z^2);
U = -mu/r + J2*R^2*mu*(3*(z/r)^2-1)/(2*r^3) + J3*mu*R^3*(5*(z/r)^3-3*(z/r))/(2*r^4);
F_x = simplify(-diff(U,x));
F_y = simplify(-diff(U,y));
F_z = simplify(-diff(U,z));

clear all 

% part b
Re = 6378;
J2 = 1082.6e-06;
J3 = 2.532153e-06;
mu = 398600.4405;

% initial conditions
a(1) = Re + 1000;
e(1) = 0.1;
i(1) = 30 *pi/180;
Omega = 270 *pi/180;
omega = 120 *pi/180;
theta = 0 *pi/180; %true anomaly

% compute initial states
p = a*(1-e^2); % orbital parameter
r = p / (1+e*cos(theta)); % distance
h = sqrt(p*mu); % angular momentum
r_P = [r*cos(theta); r*sin(theta); 0]; % position in perifocal frame
v_P = [-mu*sin(theta)/h; mu*(e+cos(theta))/h; 0];
C_G2P = [ cos(omega) sin(omega)  0;
         -sin(omega) cos(omega)  0;
          0          0           1]*...
        [ 1          0           0;
          0          cos(i)      sin(i);
          0         -sin(i)      cos(i)]*...
        [ cos(Omega) sin(Omega)  0;
         -sin(Omega) cos(Omega)  0;
          0          0           1];
C_P2G = C_G2P';
r_ECI = C_P2G*r_P;
v_ECI = C_P2G*v_P;
X0 = [r_ECI; v_ECI];

% run ode
tspan    = 0:600:24*3600*180;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
[Tp, Xp] = ode45(@(t,X) newtonJ2J3(t,X,Re,mu,J2,J3),tspan,X0, options);
timeVec = Tp;

% post processing
I = [1 0 0];J = [0 1 0];K = [0 0 1]; E = [];
for k=1:1:length(timeVec)
    r = Xp(k,1:3)';
    z = r(3);
    v = Xp(k,4:6)';
    hVec = (cross(r,v)); % Magn. of angular momentum vector per unit mass
    eVec = (1/mu)*(cross(v,hVec)-mu*r/norm(r));
    % perifocal frame unit vectors
    i_h = hVec/norm(hVec);
    i_e = eVec/norm(eVec);
    i_p = cross(i_h,i_e);
    n_0 = cross(K,i_h)/norm(cross(K,i_h));
    Omega(k) = atan2(dot(n_0,J),dot(n_0,I)); % RAAN [rad]
    i(k) = acos(dot(K,i_h));
    omega(k) = atan2(dot(i_h,cross(n_0,i_e)),dot(n_0,i_e));
    e(k) = norm(eVec);
    E(k) = 1/2*dot(v,v) - mu/norm(r); % Total energy (E=T+U) per unit mass
    E2(k) = E(k) + J2*Re^2*mu*(3*(z/norm(r))^2-1)/(2*norm(r)^3) + J3*mu*Re^3*(5*(z/norm(r))^3-3*(z/norm(r)))/(2*norm(r)^4);
end

figure(1); % 3D Position
plot3(Xp(:,1), Xp(:,2), Xp(:,3),'g-','linewidth',3); 
xlabel('x (km)','fontsize',14);ylabel('y (km)','fontsize',14); 
zlabel('z (km)','fontsize',14); 
set(gca,'fontsize',14)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
hold off

figure(2) % eccentricity
plot(timeVec/60,e,'linewidth',3)
ylabel('Eccentricity','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(3) % inclination [deg]
plot(timeVec/60,i*180/pi,'linewidth',3);
ylabel('Inclination [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(4) % RAAN [deg]
plot(timeVec/60,wrapTo2Pi(Omega)*180/pi,'linewidth',3);
ylabel('RAAN [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(5) % argument of perigee [deg]
plot(timeVec/60,wrapTo2Pi(omega)*180/pi,'linewidth',3);
ylabel('Argument of Perigee [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(6) % total energy
plot(timeVec/60,E2,'linewidth',3);
ylabel('Energy per unit mass km^2/sec^2','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(7) % nominal potential only and without perturbations
plot(timeVec/60,E,'linewidth',3);
ylabel('Energy per unit mass w/o perturbations km^2/sec^2','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(8) % nominal potential only and without perturbations
plot(timeVec/60,E,'linewidth',3);
hold on
plot(timeVec/60,E2,'linewidth',3);
ylabel('Energy km^2/sec^2','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)
legend('Energy per unit mass w/o perturbations','Energy per unit mass');

function out = newtonJ2J3(t,X,Re,mu,J2,J3)
R  = Re;
r  = norm(X(1:3));          % position vector (km)
x  = X(1);
y  = X(2);
z  = X(3);
v  = X(4:6);          % velocity vector (km/sec) 

F_x = -(2*mu*x*(x^2 + y^2 + z^2)^3 - 35*J3*R^3*mu*x*z^3 + 3*J2*R^2*mu*x*(x^2 + y^2 + z^2)^2 - 15*J2*R^2*mu*x*z^2*(x^2 + y^2 + z^2) + 15*J3*R^3*mu*x*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));
F_y = -(2*mu*y*(x^2 + y^2 + z^2)^3 - 35*J3*R^3*mu*y*z^3 + 3*J2*R^2*mu*y*(x^2 + y^2 + z^2)^2 - 15*J2*R^2*mu*y*z^2*(x^2 + y^2 + z^2) + 15*J3*R^3*mu*y*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));
F_z = -(2*mu*z*(x^2 + y^2 + z^2)^3 - 3*J3*R^3*mu*(x^2 + y^2 + z^2)^2 - 35*J3*R^3*mu*z^4 + 9*J2*R^2*mu*z*(x^2 + y^2 + z^2)^2 - 15*J2*R^2*mu*z^3*(x^2 + y^2 + z^2) + 30*J3*R^3*mu*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

Xdot(1:3)  = v;
Xdot(4:6)  = [F_x;F_y;F_z]; 
out       = Xdot(:); % convert to vector column format
end 
