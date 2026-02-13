%% HW2 PROBLEM 1
clc; clear all; close all

% create data
dt = 0.01;
timeVec=0:dt:30;
X0 = [0.5; 1];
mu = 2;

[tOut,yOut]=ode45(@(t,X) vdp_rhs(t,X,mu),timeVec,X0);

% create model
KMatrix=[];
CMatrix=[];
l=1;
dataArray = yOut(:,1);
featuresMatrixX =[];
featuresMatrixY =[];

% create feature matrices with history, one step back
for ii=2:length(dataArray)-1

featuresArrayX = [hermitePoly(dataArray(ii,1)); hermitePoly(dataArray(ii-1,1))]; % x_{n} x_{n-1} 
featuresArrayY = [hermitePoly(dataArray(ii+1,1)); hermitePoly(dataArray(ii,1))]; % x_{n+1} x_{n}
featuresMatrixX = [featuresMatrixX featuresArrayX];
featuresMatrixY = [featuresMatrixY featuresArrayY];

end

KMatrix = featuresMatrixY*pinv(featuresMatrixX);
CMatrix = dataArray(2:end-1)'*pinv(featuresMatrixX);

% % --- change initial conditions if desired
% X0 = [3; 0.1];
% % create dynamics
% [tOut,yOut]=ode45(@(t,X) vdp_rhs(t,X,mu),timeVec,X0);
% % ----

% get x1 dynamics from model (from C K X0)
x(1)=X0(1,1); % x->x1 at t0
x(2) = x(1)+dt*X0(2,1); % x->x1 at t1
for ii=2:length(dataArray)-1
x(ii+1) = CMatrix*KMatrix*[hermitePoly(x(ii)); hermitePoly(x(ii-1))];
end

figure(1)
% scatter(1:length(yOut(:,1)),yOut(:,1),'bo',LineWidth=0.1);
plot(yOut(:,1),LineWidth=1,Color='b');
hold on 
plot(x,LineWidth=1,Color='m')
xlim([1 length(yOut(:,1))])
grid on
legend('exact','prediction')
xlabel('time steps')
ylabel('x^1')
title('Van der pol Oscillator')

function out = hermitePoly(x)
out       = [2*x; 
             4*x^2-2;
             8*x^3-12*x;
             16*x^4-48*x^2+12];
end 


function out = vdp_rhs(t,X,mu)
out       = [X(2); mu*(1-X(1)^2)*X(2)-X(1)];
end 