%HW7Q1
clear all; close all; clc;

% constants
Re = 6378;
mu = 398600.4405;
J2 = 1082.6e-06;
J3 = 2.532153e-06*0;

% initial conditions
a(1) = Re + 400;
e(1) = 0.02;
i(1) = 35 *pi/180;
Omega = 90 *pi/180;
omega = 30 *pi/180;
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
tspan    = 0:600:345600;
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
    a(k) = norm(hVec)^2/mu/(1-e(k)^2);
end


figure(1) % semimajor
plot(timeVec/60,a,'linewidth',3)
ylabel('Semi-major Axis','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

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