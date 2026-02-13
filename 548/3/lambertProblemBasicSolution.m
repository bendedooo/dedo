%% Problem 5
function main
close all
clear all
departureDelay = 10 ; % [days]
mu = 1.327124e11; % Sun Gravitational parameter (km^3/s^2) 
Re  = 6378;
% Departure data
R0   = 1e8*[-1.132435294297716  -0.995058041957731   0.000049491722571]'; % Earth May 2, 2018
VE0  = [19.177994958478273 -22.489431572900006   0.000732470689783]';     % Earth velocity at departure
JD0  = 2.4582405e6; % Julian date of departure

if departureDelay > 0
X0 = [R0; VE0];
tspan    = 0:600:24*3600*departureDelay;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
[Tp, Xp] = ode45(@(t,X) satdyn(t,X),tspan,X0, options);
R0 = Xp(end,1:3)';
VE0 = Xp(end,4:6)';
JD0 = JD0 + departureDelay;
end

% Arrival data
Rf1  = 1e8*[1.433197554726618   1.669283013110051  -0.000214501948468]';  % Mars Jan 15, 2019
VM1  = [-17.461633983252142  17.849763278562293   0.802958937672343]'; 
JD1  = 2.4584985e6; % Julian date of arrival

Ri     = R0;
Rf     = Rf1;
JDf    = JD1; 
VEarth = VE0(:);
VMars  = VM1(:);

TOFdays = JD1 - JD0;

TOF = (TOFdays)*3600*24; % convert days into TOF in sec

[p, Vi, Vf, f, g, fd, gd] = lambert(mu, Ri, Rf, TOF);

totalDeltaV = norm(Vi(:) - VEarth(:)) + norm(Vf(:) - VMars(:))


X0 = [Ri; Vi];
tspan    = 0:600:TOF;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
[Tp, Xp] = ode45(@(t,X) satdyn(t,X),tspan,X0, options);
timeVec = Tp;

figure(1); % 3D Position
plot3(Xp(:,1), Xp(:,2), Xp(:,3),'g-','linewidth',3); 
xlabel('x (km)','fontsize',16);ylabel('y (km)','fontsize',16); 
zlabel('z (km)','fontsize',16); 
hold on;



X0 = [Ri; VEarth];
tspan    = 0:600:TOF*2;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
[TpE, XpE] = ode45(@(t,X) satdyn(t,X),tspan,X0, options);
timeVec = TpE;
plot3(XpE(:,1), XpE(:,2), XpE(:,3),'b-','linewidth',3); 
hold on

X0 = [Rf; VMars];
tspan    = 0:600:TOF*4;
options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
[TpM, XpM] = ode45(@(t,X) satdyn(t,X),tspan,X0, options);
timeVec = TpM;
plot3(XpM(:,1), XpM(:,2), XpM(:,3),'r-','linewidth',3); 

% set(gca,'fontsize',14)
% [XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
% surf(XS*Re, YS*Re, ZS*Re);
% axis equal

scatter3(XpE(1,1), XpE(1,2), XpE(1,3),'b*','linewidth',3); 
scatter3(XpM(1,1), XpM(1,2), XpM(1,3),'r*','linewidth',3); 

legend('Transfer Orbit', 'Earth Orbit', 'Mars Orbit', 'Earth Initial Position', 'Mars Initial Position');
% hold off
return






% Brute force solver for Lambert problem - written by ilya@umich.edu
function [p, Vi, Vf, f, g, fd, gd] = lambert(mu, Ri, Rf, TOF)
global A muG t ri rf
muG = mu;
t  = TOF;
ri  = norm(Ri); rf = norm(Rf);
q   = dot([0;0;1],cross(Ri,Rf));
phi = acos(dot(Ri,Rf)/(ri*rf)); 
if (q < 0),
    phi=2*pi - phi;
end
A = sqrt(ri*rf/(1-cos(phi)))*sin(phi);
z = bisection(@Fnon, 0.001, 20);
[res, x, y] = Fnon(z);
p = ri*rf*(1 - cos(phi))/y;
f = 1 - y/ri;
g = A*sqrt(y/mu);
gd = 1      - y/rf;
fd = (f*gd  - 1)/g;
Vi = (Rf    - f*Ri)/g;
Vf = (gd*Rf - Ri)/g; 
return;

function [res,x,y] = Fnon(z)
global A muG t ri rf 
[C, S] = stumpff(z,20);
c3 = S; c2 = 0; c1 = A*sqrt(C); c0 = -sqrt(muG)*t;
roots_all  = roots([c3,c2,c1,c0]);
x = roots_all(find(abs(imag(roots_all))<1e-6));
y = C*x^2;
zout = (1- ( sqrt(C)*(ri + rf - y)/A  ))/S;
res = z - zout;
return;

function [C,S]=stumpff(z, n)
C = 0.5;
S = 1/6;
for i=1:1:n,
    C=C+(-1)^i*z^i/factorial(2*(i+1));
    S=S+(-1)^i*z^i/factorial(2*i+3);
end;

return

function p = bisection(f,a,b)
if f(a)*f(b)>0
    disp('No Solution Found: Pick different interval [a,b]')
else
    p = (a + b)/2;
    err = abs(f(p));
    while err > 1e-12
        if f(a)*f(p)<0
            b = p;
        else
            a = p;
        end
        p = (a + b)/2;
        err = abs(f(p));
    end
end
return

function Xdot = satdyn(t,X)
 
mu = 1.327124e11;        % gravitational constant for Earth (km^3/s^2) 
r  = X(1:3);          % position vector (km)
v  = X(4:6);          % velocity vector (km/sec) 
rn = norm(r);
Xdot(1:3)  = v;
Xdot(4:6)  = -mu/rn^3*r ; 
Xdot       = Xdot(:); % convert to vector column format
 
return;