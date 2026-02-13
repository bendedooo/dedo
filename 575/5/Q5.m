%% part a
clear all; clc; close all;
tf = 1;
tspan = [tf:-0.01:0]; % Note time is running backward 
x0 = [1; 0];
% define the system

R = 1;
Q = zeros(2);
A = [0 1; 0 0];
B = [0; 1];
muVec = [0.01 1 10 100];
for i=1:length(muVec)
mu = muVec(i);
Sf = eye(2).*mu;

% Initialize P(t1)
Yt1 = [Sf(1,1); Sf(1,2); Sf(2,1); Sf(2,2)]; 
tspan = [tf:-0.01:0]; % Note time is running backward 

[T,Y] = ode45(@(t,Y) ric(t,Y,A, B, Q, R), tspan, Yt1 ); % Integrate backward-in-time Riccati equation
inds = length(T):-1:1; % reverse time

Tout = T(inds);
Yout = Y(inds,:);
% Simulate closed-loop system
tspan = [0:0.01:1]; 
[T, X] = ode45(@(t,x) clp(t, x,A, B, Q, R, Tout, Yout), tspan, x0 );
U=[];
for m=1:1:length(T)
    [xdot,u] = clp(T(m),X(m,:)',A, B, Q, R, Tout, Yout);
    U = [U;u];
end
figure(1)
subplot(4,1,i)
plot(T,X(:,1),T,X(:,2),T,U,'linewidth',2);
if i == 1
legend('x_1','x_2','u','Location','best')
end
xlabel('Time [sec]')
ylabel('x,u'); set(gca,'fontsize',12)
title("Problem 5A: Optimal Control and State Trajectories with \mu =" + mu)
set(gcf,'units','points','position',[10,10,1000,1000])
% X(end,:) % states change the most as mu increases, penalty in final state is larger in the cost function,  requires bigger change to minimize the cost
% % penalty on missing the desired state is high, so it tries bigger control? % mu increases  increases the weight on the last state in the cost function, u gets bigger to minimize J
end
%% part b

% Initialize P(t1)
% Simulate closed-loop system
tspan = [0:0.01:tf]; 
[T, X] = ode45(@(t,x) clpb(t, x,A, B), tspan, x0 );
U=[];
for m=1:1:length(T)
    [xdot,u] = clpb(T(m),X(m,:)',A, B);
    U = [U;u];
end
figure(i+1)
plot(T,X(:,1),T,X(:,2),T,U,'linewidth',2);
legend('x_1','x_2','u')
xlabel('Time [sec]')
ylabel('x,u'); set(gca,'fontsize',12)
title("Problem 5B: Optimal Control and State Trajectories Hard Terminal LQ")

%hard lq requires more control effort, same as soft as mu increases 
%% functions

function Ydot = ric(t,Y,A, B, Q, R)
P = [Y(1),Y(2);Y(3),Y(4)]; % re-compose P matrix from vector
Pdot = -A'*P - P*A + P*B*inv(R)*B'*P - Q; % compute P-dot
Ydot = [Pdot(1,1);Pdot(1,2);Pdot(2,1);Pdot(2,2)]; % organize into a vector
end

function [xdot, u] = clp(t, x,A, B, Q, R, Tout, Yout)
p11 = interp1(Tout, Yout(:,1), t); % 1D interpolation 
p12 = interp1(Tout, Yout(:,2), t);
p21 = interp1(Tout, Yout(:,3), t);
p22 = interp1(Tout, Yout(:,4), t);
P = [p11, p12; p21, p22];
u = -inv(R)*B'*P*x;
xdot = A*x + B*u;
end

function [xdot, u] = clpb(t, x,A, B)
u = -6 + 12*t;
xdot = A*x + B*u;
end

