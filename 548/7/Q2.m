%HW7Q2
clear all; close all; clc;

% constants
mu  = 398600.4405;

% initial conditions
a0 = 7000;
e0 = 0.1;
i0 = pi/6;
Omega0 = pi;
omega0 = 0;
theta0 = 0; %true anomaly

%desired conditions
ades = 10000;
edes = 0.2;

tspan    = 0:600:40*60*60;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
X0 = [a0; e0; i0; Omega0; omega0; theta0]; 
[Tp, Xp] = ode45(@(t,X) GVE(t,X,mu, ades, edes),tspan,X0, options);
timeVec = Tp;

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
[oout, S(i), T(i), W(i)] = GVE(tspan(i),Xp(i,:),mu,ades, edes);
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

function [out, Sout, Tout, Wout] = GVE(t,X,mu,ades, edes)

a     = X(1);
e     = X(2);
i     = X(3);
Omega = X(4);
omega = X(5);
theta = X(6); %true anomaly

p     = a*(1-e^2);
r     = p/(1+e*cos(theta));
cpsi  = (1/e)*(1-r/a);

gamma1 = 1e-05;
gamma2 = 1e03;

% control law
K = [1e-06  0; 0 1e-06];
B = [ (2*a^2/sqrt(mu*p))*(e*sin(theta))   (2*a^2/sqrt(mu*p))*(p/r);
      (p*sin(theta)/sqrt(mu*p))             (p*(cpsi+cos(theta))/sqrt(mu*p))];
F = -K*B'*[gamma1*(a-ades); gamma2*(e-edes)];
S = F(1);
T = F(2);
W = 0;

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
