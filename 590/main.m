% sorular
% - benim update rate'im ne olacak?, odeyi onunla calistiracagim
% - how did we choose the controller parameters, loveras paper??
% - II kismi neden uymuyoo ya farkli orbitlerde
% - check the figures 17 abcd with control
% - save fig as pdf yaparken a4 olmasin
% - miguelin dedigine bak



clear all; close all; clc;
tic
%% options/flags
% for plots
FlagOrbitPlot       = 1;
FlagAngVelPlot      = 1;
FlagMagFieldPlot    = 0;
FlagAxisPlot        = 0;
FlagEulerPlot       = 0;
flagQuaternionPlot  = 1;
flagPlotDipMoment   = 1;

% options
FlagMagField        = 1; 
FlagControl         = 2; % 0: no control, 1: saturated, 2: unsaturated

orbitRep            = 20;
constRot            = 0;
%% constants
mu  = 398600.4405;
Re  = 6378.2;

%% initial conditions
a0 = Re+450;
e0 = 0;
i0 = 0;
Omega0 = 0*pi/180;
omega0 = 0*pi/180;
theta0 = 0; %true anomaly

epoch = '1-Jan-2016 00:00:00';% ?? I.c icin tam olarak hangi zamani kullaniyoruz [1, 1, 2016, 12, 00, 00]; % [DD MM YYYY Hours Min Sec]
% epochJ2000 = '1-Jan-2000 12:00:00';
TimeStamp0 = datetime(epoch);
period     = 2*pi*sqrt(a0^3/mu);
tspan      = 0:0.5:period*orbitRep;

% create the orbit
% compute initial states
p = a0*(1-e0^2); % orbital parameter
r = p / (1+e0*cos(theta0)); % distance
h = sqrt(p*mu); % angular momentum
r_P = [r*cos(theta0); r*sin(theta0); 0]; % position in perifocal frame
v_P = [-mu*sin(theta0)/h; mu*(e0+cos(theta0))/h; 0];
C_I2P = Angle2DCM_313(Omega0,i0,omega0);
C_P2I = C_I2P';
r_ECI_0 = C_P2I*r_P;
v_ECI_0 = C_P2I*v_P;

%% adjust
omega_B_I_B_0 = [0.01 0 0]'; %rad/sec
% q_0 = [sin(-pi/2/2).*[0 1 0] cos(-pi/2/2)]';% [q_v q_4][cos(-pi/2/2) sin(-pi/2/2).*[0 1 0]]';
q_0 = [0 0 1 0]'; % [q_v q_4]'

X0 = [r_ECI_0; v_ECI_0; omega_B_I_B_0; q_0];

%% run ode
% options  = odeset('AbsTol',1e-8,'RelTol',1e-9);
options  = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
[Torb, Xorb] = ode113(@(t,X) scDyn13(t,X,TimeStamp0,FlagMagField,FlagControl,constRot),tspan,X0, options);
timeVec = Torb;

%% plots
plots
toc