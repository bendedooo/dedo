clear all; clc; close all;
% find a1 and a2
tf=100;
aVec = [0.1 11]';
fun = @myFun;
options = optimoptions('fsolve');%,'Display','iter');
aVec = fsolve(fun,aVec,options);
a1 = aVec(1)
a2 = aVec(2)
p_40 = a2+a1*tf
p_20 = a1

% run ODE with a1 and a2
tf=100;
dt=0.1;
timeVec=0:dt:tf;
states0 = [0 0 0 0]';
odeFunc = @(t,x) odefcn(t,x,aVec);
optPos = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);
for i = 1:length(tOut)
    ti = tOut(i);
    uOut(i) = atan2(a1*ti-a1*tf-a2,1);
end

%% plots
close all

figure(1) % u vs t
plot(tOut,uOut*180/pi,'LineWidth',2)
xlabel('t'); ylabel('u^*'); grid on;
title('Problem 2B Optimal Control u^*')


figure(2) % v1 vs t
plot(tOut,statesOut(:,3),'LineWidth',2)
xlabel('t'); ylabel('v^*_1'); grid on;
title('Problem 2B v^*_1')


figure(3) % y vs x
plot(statesOut(:,1),statesOut(:,2),'LineWidth',2)
hold on
scatter(statesOut(1,1),statesOut(2,1),'filled','r')
xlabel('x'); ylabel('y'); grid on;
title('Problem 2B Trajectory')
legend('Optimal Rocket Trajectory','Rocket Initial Point')


%% functions
function dydt = odefcn(t,x,aVec)
tf=100;
% dt = 0.1;
% t = k*dt;
a1 = aVec(1);
a2 = aVec(2);
u = atan2(a1*t-a1*tf-a2,1);
dydt = zeros(2,1);
a = 100;
dydt(1) = x(3);
dydt(2) = x(4);
dydt(3) = a*cos(u);
dydt(4) = a*sin(u);
end

function F = myFun(aVec)
h = 8000;
tf=100;
dt=0.1;
timeVec=0:dt:tf;

states0 = [0 0 0 0]';

odeFunc = @(t,x) odefcn(t,x,aVec);
optPos = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);

% u = atan2(aVec(1)*(tf-dt)-aVec(1)*tf-aVec(2),1);
%
% v2 = statesOut(end-1,4) + a*sin(u)*dt;
% y  = statesOut(end-1,2) + (statesOut(end-1,4) + a*sin(u)*dt)*dt;

v2 = statesOut(end,4);
y  = statesOut(end,2);

F(1) = y - h;
F(2) = v2;
end

