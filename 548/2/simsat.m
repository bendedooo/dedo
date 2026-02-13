function simsat
 
Re      = 6378;      % Earth Radius (km)
mu      = 398600.4405;  % Gravitational constant for Earth (km^3/s^2) 
 
% Initial position vector in ECI frame (units: km)
R0 = [6978; 0; 0]; % [-8.058685464116713e+002;-4.135050989428283e+003;5.071447347723579e+003];  
 
% Initial velocity vector in ECI frame (units: km/sec)
V0 = [0; 5.8787; 5.8787]; % [6.04068079277917;3.52422658592595;3.83339310219652];     
 
% Initial state [position vector on top of velocity vector] 
X0 = [R0; V0]; 
 
t0     = 0;          % initial time, sec
tstep  = 20;         % time step to sample the output, sec
tfinal = 108000 ; % 96.2*60*5;  % final time, sec, approximately 5 orbits
 
% Integrate ODEs using Runge-Kutta 45 method
tspan    = t0:tstep:tfinal;
options = odeset('AbsTol',1e-6,'RelTol',1e-8); 
[Tp, Xp] = ode45(@satdyn,tspan,X0,options);
timeVec = Tp;XTraj = Xp;

figure(1);
plot3(Xp(:,1), Xp(:,2), Xp(:,3),'g-','linewidth',3); 
xlabel('x (km)','fontsize',16);ylabel('y (km)','fontsize',16); 
zlabel('z (km)','fontsize',16); 
set(gca,'fontsize',14)
 
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
hold on;
surf(XS*Re, YS*Re, ZS*Re);
% axis([-1e4, 1e4, -1e4, 1e4, -1e4, 1e4])
axis equal
hold off
 
% plot altitude above the Earth, km
figure(2);
plot(Tp/60, sqrt(Xp(:,1).^2 + Xp(:,2).^2 + Xp(:,3).^2) - Re, 'linewidth',3);
ylabel('Altitude r - R_e (km)','fontsize',16)
xlabel('t (min)','fontsize',16)
set(gca,'fontsize',16)

E = []; h=[]; e_h=[];
for k=1:1:length(timeVec)
    r = XTraj(k,1:3)';
    v = XTraj(k,4:6)';
    E(k) = 1/2*dot(v,v) - mu/norm(r); % Total energy (E=T+U) per unit mass
    hVec = (cross(r,v)); % Magn. of angular momentum vector per unit mass
    h(k) = norm(hVec); % Magn. of angular momentum vector per unit mass
    e_h(:,k)        = hVec/h(k);
end

figure(5);
plot(timeVec/60,E,'linewidth',3);
ylabel('Energy per unit mass km^2/sec^2','fontsize',16) 
xlabel('t (min)','fontsize',16)
set(gca,'fontsize',16)

figure(3);
plot(timeVec/60,h,'linewidth',3);
ylabel('Magn. of ang. mom. per unit mass km^2/sec','fontsize',14);
xlabel('t (min)','fontsize',16)
set(gca,'fontsize',16)

figure(4)
subplot(3,1,1)
plot(timeVec/60,e_h(1,:),'linewidth',3);
ylabel('e_h_x','fontsize',14);
xlabel('t (min)','fontsize',16)

subplot(3,1,2)
plot(timeVec/60,e_h(2,:),'linewidth',3);
ylabel('e_h_y','fontsize',14);
xlabel('t (min)','fontsize',16)

subplot(3,1,3)
plot(timeVec/60,e_h(3,:),'linewidth',3);
ylabel('e_h_z','fontsize',14);
xlabel('t (min)','fontsize',16)
sgtitle('Components of e_h') 
return

function Xdot = satdyn(t,X)
 
mu = 398600.4405;        % gravitational constant for Earth (km^3/s^2) 

part = 'c';

if part =='a'
     P =[0;0;0];
elseif part =='b'
     P =[5e-4;0;0];
elseif part =='c'
     P =[0;2e-5;0];
elseif part =='d'
     P =[0;0;2e-3];    
end

r  = X(1:3);          % position vector (km)
v  = X(4:6);          % velocity vector (km/sec) 
h  = cross(r,v);

rn         = norm(r); % distance (km)
hn         = norm(h);
e_r        = r/rn;
e_h        = h/hn;
e_theta    = cross(e_h,e_r);


Xdot(1:3)  = v;
Xdot(4:6)  = -mu/rn^3*r + P(1)*e_r + P(2)*e_h + P(3)*e_theta; 
 


Xdot       = Xdot(:); % convert to vector column format
 
return;