%HW7Q3
clear all; close all; clc;

% constants
mu  = 398600.4405;
Re  = 6378.2;
J2  = 1082.64e-06;
AJ2 = 0.5*J2*Re^2;

% initial conditions
a0 = 7000;
e0 = 0.1;
i0 = pi/6;
Omega0 = pi;
omega0 = 0;
theta0 = 0; %true anomaly

tspan    = 0:600:40*60*60;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
X0 = [a0; e0; i0; Omega0; omega0; theta0];
[Tp, Xp] = ode45(@(t,X) GVE(t,X,mu),tspan,X0, options);
timeVec = Tp;

% post processing
for k=1:1:length(timeVec)
    a     = Xp(k,1);
    e     = Xp(k,2);
    i     = Xp(k,3);
    Omega = Xp(k,4);
    omega = Xp(k,5);
    theta = Xp(k,6); %true anomaly
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
    rOut(:,k) = r_ECI;

end
figure(1) % semimajor
plot(timeVec/60,Xp(:,1),'linewidth',3)
ylabel('Semi-major Axis','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(2) % eccentricity
plot(timeVec/60,Xp(:,2),'linewidth',3)
ylabel('Eccentricity','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

% figure(8) % inclination [deg]
% plot(timeVec/60,Xp(:,3)*180/pi,'linewidth',3);
% ylabel('Inclination [deg]','fontsize',14)
% xlabel('t (min)','fontsize',14)
% set(gca,'fontsize',14)
% figure(9) % RAAN [deg]
% plot(timeVec/60,wrapTo2Pi(Xp(:,4))*180/pi,'linewidth',3);
% ylabel('RAAN [deg]','fontsize',14)
% xlabel('t (min)','fontsize',14)
% set(gca,'fontsize',14)

figure(3) % argument of perigee [deg]
plot(timeVec/60,wrapTo2Pi(Xp(:,5))*180/pi,'linewidth',3);
ylabel('Argument of Perigee [deg]','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)
for i=1:length(tspan)
    [oout, S(i), T(i), W(i)] = GVE(tspan(i),Xp(i,:),mu);
end


figure(4) %
plot(timeVec/60,S,'linewidth',3)
ylabel('S','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(5) %
plot(timeVec/60,T,'linewidth',3)
ylabel('T','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(6) %
plot(timeVec/60,W,'linewidth',3)
ylabel('W','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)


figure(7); % 3D Position
plot3(rOut(1,:), rOut(2,:), rOut(3,:),'g-','linewidth',3);
xlabel('x (km)','fontsize',14);ylabel('y (km)','fontsize',14);
zlabel('z (km)','fontsize',14);
set(gca,'fontsize',14)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
title('3D Trajectory GVE')
hold off

%% part b
% initial conditions
a0 = 7000;
e0 = 0.1;
i0 = pi/6;
Omega0 = pi;
omega0 = 0;
theta0 = 0; %true anomaly

% compute initial states
p = a0*(1-e0^2); % orbital parameter
r = p / (1+e0*cos(theta0)); % distance
h = sqrt(p*mu); % angular momentum
r_P = [r*cos(theta0); r*sin(theta0); 0]; % position in perifocal frame
v_P = [-mu*sin(theta0)/h; mu*(e0+cos(theta0))/h; 0];
C_G2P = [ cos(omega0) sin(omega0)  0;
    -sin(omega0) cos(omega0)  0;
    0          0           1]*...
    [ 1          0           0;
    0          cos(i0)      sin(i0);
    0         -sin(i0)      cos(i0)]*...
    [ cos(Omega0) sin(Omega0)  0;
    -sin(Omega0) cos(Omega0)  0;
    0          0           1];
C_P2G = C_G2P';
r_ECI = C_P2G*r_P;
v_ECI = C_P2G*v_P;
X0 = [r_ECI; v_ECI];

% run ode
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
figure(8) % semimajor
plot(timeVec/60,a,'linewidth',3)
ylabel('Semi-major Axis','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(9) % eccentricity
plot(timeVec/60,e,'linewidth',3)
ylabel('Eccentricity','fontsize',14) 
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(10) % argument of perigee [deg]
plot(timeVec/60,(omega)*180/pi,'linewidth',3);
ylabel('Argument of Perigee [deg]','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

for i=1:length(tspan)
    [oout, S(i), T(i), W(i)] = ECIJ2(tspan(i),Xp(i,:),mu,AJ2);
end


figure(11) %
plot(timeVec/60,S,'linewidth',3)
ylabel('S','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(12) %
plot(timeVec/60,T,'linewidth',3)
ylabel('T','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(13) %
plot(timeVec/60,W,'linewidth',3)
ylabel('W','fontsize',14)
xlabel('t (min)','fontsize',14)
set(gca,'fontsize',14)

figure(14); % 3D Position
plot3(Xp(:,1), Xp(:,2), Xp(:,3),'g-','linewidth',3);
xlabel('x (km)','fontsize',14);ylabel('y (km)','fontsize',14);
zlabel('z (km)','fontsize',14);
set(gca,'fontsize',14)
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
axis equal
title('3D Trajectory Conventional EOM')
hold off


% functions
function [out, Sout, Tout, Wout] = GVE(t,X,mu)

a     = X(1);
e     = X(2);
i     = X(3);
Omega = X(4);
omega = X(5);
theta = X(6); %true anomaly

p     = a*(1-e^2);
r     = p/(1+e*cos(theta));
cpsi  = (1/e)*(1-r/a);
rdot  = sqrt(mu/p)*e*sin(theta);

% C_SMF2P = [ cos(theta) sin(theta)  0;
%     -sin(theta) cos(theta)  0;
%     0          0           1];
% C_P2SMF = C_SMF2P';
h = sqrt(p*mu);
% v_P = [-mu*sin(theta)/h; mu*(e+cos(theta))/h; 0];
% v_SMF = C_P2SMF*v_P;
v_SMF = [rdot; h/r; 0];

e_v   = v_SMF/norm(v_SMF);

F  = (1e-05)*e_v;

S = F(1);
T = F(2);
W = F(3);

Xdot(1) = (2*a^2/sqrt(mu*p))*(e*sin(theta)*S+p*T/r);
Xdot(2) = (p*sin(theta)/sqrt(mu*p)) * S + (p*(cpsi+cos(theta))/sqrt(mu*p)) * T;
Xdot(3) = (r*cos(theta+omega)/sqrt(mu*p)) * W;
Xdot(4) = (r*sin(theta+omega)/(sqrt(mu*p)*sin(i))) * W;
Xdot(5) = -(p*cos(theta)/(e*sqrt(mu*p))) * S + ((r+p)*sin(theta)/(e*sqrt(mu*p))) * T - (r*sin(theta+omega)*cot(i)/sqrt(mu*p)) * W;
Xdot(6) = sqrt(mu*p)/r^2 + (p*cos(theta)/(e*sqrt(mu*p))) * S - ((p+r)*sin(theta)/(e*sqrt(mu*p))) * T;

out     = Xdot(:); % convert to vector column format
Sout    = S;
Tout    = T;
Wout    = W;
end

function [out, Sout, Tout, Wout] = ECIJ2(t,X,mu,AJ2)
r  = norm(X(1:3));          % position vector (km)
x  = X(1);
y  = X(2);
z  = X(3);
v  = X(4:6);          % velocity vector (km/sec)

e_v   = v/norm(v);

Fpot   = -mu/r^3*X(1:3); 
Fpert  = (1e-05)*e_v;

F_x = Fpot(1)+Fpert(1);
F_y = Fpot(2)+Fpert(2);
F_z = Fpot(3)+Fpert(3);

% compute SMF unit vectors
h          = (cross(X(1:3),v));
rn         = norm(X(1:3));
hn         = norm(h);
e_r        = X(1:3)/rn;
e_h        = h/hn;
e_theta    = cross(e_h,e_r);

Sout = dot(Fpert ,e_r);
Tout = dot(Fpert ,e_theta);
Wout = dot(Fpert ,e_h);


Xdot(1:3)  = v;
Xdot(4:6)  = [F_x;F_y;F_z];

out       = Xdot(:); % convert to vector column format




end