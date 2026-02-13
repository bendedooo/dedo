%% PROBLEM 3

function simsat
 
Re      = 6378;      % Earth Radius (km)
mu      = 398600.4405;  % Gravitational constant for Earth (km^3/s^2) 
 
% Initial position vector in ECI frame (units: km)
R0 = [0; Re+800; 1000]; % [-8.058685464116713e+002;-4.135050989428283e+003;5.071447347723579e+003];  
 
% Initial velocity vector in ECI frame (units: km/sec)
V0 = [8.9423; 0; -0.3]; % [6.04068079277917;3.52422658592595;3.83339310219652];     
 
% Initial state [position vector on top of velocity vector] 
X0 = [R0; V0]; 
 
t0     = 0;          % initial time, sec
tstep  = 20;         % time step to sample the output, sec
tfinal = 540000;%00 ; % 96.2*60*5;  % final time, sec, approximately 5 orbits
 
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
    U = -(1/6)*39.86e04*( 8/(sqrt(r(1)^2+r(2)^2+r(3)^2)) - 1/(sqrt((r(1)+3189)^2+r(2)^2+r(3)^2)) - 1/(sqrt((r(1)-3189)^2+r(2)^2+r(3)^2)) );
    E(k) = 1/2*dot(v,v) +U; % Total energy (E=T+U) per unit mass
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

subplot(3,1,2)
plot(timeVec/60,e_h(2,:),'linewidth',3);

subplot(3,1,3)
plot(timeVec/60,e_h(3,:),'linewidth',3);

return

function Xdot = satdyn(t,X)
 
mu = 398600.4405;        % gravitational constant for Earth (km^3/s^2) 



r  = X(1:3);          % position vector (km)
v  = X(4:6);          % velocity vector (km/sec) 
h  = cross(r,v);

x = r(1);
y = r(2);
z = r(3);
% syms x y z
% U = (1/6)*39.86e04*( 8/(sqrt(x^2+y^2+z^2)) + 1/(sqrt((x+3189)^2+y^2+z^2)) + 1/(sqrt((x-3189)^2+y^2+z^2)) );
% 
% fx = simplify(diff(U,x));
% fy = simplify(diff(U,y));
% fz = simplify(diff(U,z));
% fx = (4565263904494933*(2*x - 6378))/(137438953472*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*x)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*(2*x + 6378))/(137438953472*((x + 3189)^2 + y^2 + z^2)^(3/2));
% fy = (4565263904494933*y)/(68719476736*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*y)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*y)/(68719476736*((x + 3189)^2 + y^2 + z^2)^(3/2));
% fz = (4565263904494933*z)/(68719476736*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*z)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*z)/(68719476736*((x + 3189)^2 + y^2 + z^2)^(3/2));

F_x = (4565263904494933*(2*x - 6378))/(137438953472*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*x)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*(2*x + 6378))/(137438953472*((x + 3189)^2 + y^2 + z^2)^(3/2));
F_y = (4565263904494933*y)/(68719476736*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*y)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*y)/(68719476736*((x + 3189)^2 + y^2 + z^2)^(3/2));
F_z = (4565263904494933*z)/(68719476736*((x - 3189)^2 + y^2 + z^2)^(3/2)) - (4565263904494933*z)/(8589934592*(x^2 + y^2 + z^2)^(3/2)) + (4565263904494933*z)/(68719476736*((x + 3189)^2 + y^2 + z^2)^(3/2));

% F_x = single(subs(fx));
% F_y = single(subs(fy));
% F_z = single(subs(fz));


Xdot(1:3)  = v;
Xdot(4:6)  = [F_x;F_y;F_z]; 


Xdot       = Xdot(:); % convert to vector column format
 
return;