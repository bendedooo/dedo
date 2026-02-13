%HW7Q1
clear all; close all; clc;

% constants
Re  = 6378.2;
mu  = 398600.4405;
J2  = 1082.64e-06;
AJ2 = 0.5*J2*Re^2;

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
[Tp, Xp] = ode45(@(t,X) ECIJ2(t,X,mu,AJ2),tspan,X0, options);
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
% plots
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




%% part b
a0 = Re + 400;
e0 = 0.02;
i0 = 35 *pi/180;
Omega0 = 90 *pi/180;
omega0 = 30 *pi/180;
theta0 = 0 *pi/180; %true anomaly

X0 = [a0; e0; i0; Omega0; omega0; theta0]; 
[Tp, Xp] = ode45(@(t,X) GVE(t,X,mu,Re, J2),tspan,X0, options);
timeVec = Tp;


figure(6) % semimajor
plot(timeVec/60,Xp(:,1),'linewidth',3)
ylabel('Semi-major Axis','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(7) % eccentricity
plot(timeVec/60,Xp(:,2),'linewidth',3)
ylabel('Eccentricity','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(8) % inclination [deg]
plot(timeVec/60,Xp(:,3)*180/pi,'linewidth',3);
ylabel('Inclination [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)


figure(9) % RAAN [deg]
plot(timeVec/60,wrapTo2Pi(Xp(:,4))*180/pi,'linewidth',3);
ylabel('RAAN [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(10) % argument of perigee [deg]
plot(timeVec/60,wrapTo2Pi(Xp(:,5))*180/pi,'linewidth',3);
ylabel('Argument of Perigee [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

%% part c
n = sqrt(mu/a0^3);

omegadot = -3*J2*n*(1-5*cos(i0)^2)*Re^2/(4*a0^2*(1-e0^2)^2);
Omegadot = -3*J2*n*cos(i0)*Re^2/(2*a0^2*(1-e0^2)^2);

omega_avg = (omega0 + omegadot*tspan(end))*180/pi
Omega_avg = (Omega0 + Omegadot*tspan(end))*180/pi

figure(9) % RAAN [deg]
plot(timeVec/60,wrapTo2Pi(Xp(:,4))*180/pi,'linewidth',3);
hold on
plot(timeVec/60,wrapTo2Pi(Omega0+Omegadot*tspan)*180/pi,'linewidth',3);
ylabel('RAAN [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)
legend('ODE','Averaged')

figure(10) % argument of perigee [deg]
plot(timeVec/60,wrapTo2Pi(Xp(:,5))*180/pi,'linewidth',3);
hold on
plot(timeVec/60,wrapTo2Pi(omega0+omegadot*tspan)*180/pi,'linewidth',3);
ylabel('Argument of Perigee [deg]','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)
legend('ODE','Averaged')



%% functions
function out = ECIJ2(t,X,mu,AJ2)
r  = norm(X(1:3));          % position vector (km)
x  = X(1);
y  = X(2);
z  = X(3);
v  = X(4:6);          % velocity vector (km/sec) 

F_x = mu*(-x/r^3+AJ2*(15*x*z^2/r^7-3*x/r^5));
F_y = mu*(-y/r^3+AJ2*(15*y*z^2/r^7-3*y/r^5));
F_z = mu*(-z/r^3+AJ2*(15*z^3/r^7  -9*z/r^5));

Xdot(1:3)  = v;
Xdot(4:6)  = [F_x;F_y;F_z]; 
out       = Xdot(:); % convert to vector column format
end 


function out = GVE(t,X,mu,Re,J2)

a     = X(1);
e     = X(2);
i     = X(3);
Omega = X(4);
omega = X(5);
theta = X(6); %true anomaly

p     = a*(1-e^2);
r     = p/(1+e*cos(theta));
cpsi  = (1/e)*(1-r/a);

% % compute SMF unit vectors
% h          = (cross(r,v)); 
% rn         = norm(r);
% hn         = norm(h);
% e_r        = r/rn;
% e_h        = h/hn;
% e_theta    = cross(e_h,e_r);

coeff = -3*mu*J2*Re^2/(2*r^4);
S     = coeff * (1-3*sin(i)^2*sin(theta+omega)^2);
T     = coeff * 2*sin(i)^2*sin(theta+omega)*cos(theta+omega);
W     = coeff * 2*sin(i)*sin(theta+omega)*cos(i);

Xdot(1) = (2*a^2/sqrt(mu*p))*(e*sin(theta)*S+p*T/r);
Xdot(2) = (p*sin(theta)/sqrt(mu*p)) * S + (p*(cpsi+cos(theta))/sqrt(mu*p)) * T;
Xdot(3) = (r*cos(theta+omega)/sqrt(mu*p)) * W;
Xdot(4) = (r*sin(theta+omega)/(sqrt(mu*p)*sin(i))) * W;
Xdot(5) = -(p*cos(theta)/(e*sqrt(mu*p))) * S + ((r+p)*sin(theta)/(e*sqrt(mu*p))) * T - (r*sin(theta+omega)*cot(i)/sqrt(mu*p)) * W;
Xdot(6) = sqrt(mu*p)/r^2 + (p*cos(theta)/(e*sqrt(mu*p))) * S - ((p+r)*sin(theta)/(e*sqrt(mu*p))) * T;

out     = Xdot(:); % convert to vector column format
end
