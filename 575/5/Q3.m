
%% part a
clear all; clc; close all;
A = [-0.0257, 0.0130, -0.322; 1.260, -1.765, 0; 0, 1, 0];
B = [0.0860; -7.41; 0];
C = [1, 0, 0];

% check controllability and observability
isControllable = rank([B, A*B, A*A*B]) == length(A) % = rank(ctrb(A,B))
isObsevable    = rank([C; C*A; C*A*A]) == length(A) % = rank(obsv(A,C))

%% part b
qVec = [0, 1, 100, 1000, 10000];
R = 1;
for i = 1:length(qVec)
    q = qVec(i);
    Q = q*C'*C;
    [K,~,~] = lqr(A,B,Q,R); % u = -kx
    KMatrix(i,:) = -K;
    clSysA = A+B*KMatrix(i,:);
    eValclSys(:,i) = eig(clSysA);
    isStableclSys(i) = prod(real(eValclSys(:,i))<0);% (max(vecnorm(eValclSys(:,i),2,2))<=1);
end
%% part c
K = KMatrix(2,:);
tspan = [0:0.01:10];
x0 = [1;0;0];
[T,X]=ode45( @(t,x) myclp(t,x, A, B, K), tspan, x0 );
U=[];
for m=1:1:length(T)
    [xdot,u] = myclp(T(m),X(m,:)',A,B,K);
    U = [U;u];
end
figure(1);
plot(T,X(:,1),T,X(:,2),T,X(:,3),T,U,'linewidth',2);
legend('x_1','x_2','x_3','u')
xlabel('Time [sec]')
ylabel('x,u'); set(gca,'fontsize',12)
title('Problem 3C: Closed-Loop Time Responses')
%% part d
global r;
r = 1;
x0 = [0,0,0]';
R = 1;
Q1 = 100;
Q = [Q1 zeros(1,3); zeros(3,1) zeros(3,3)];
D = 0; 
At = [0 C; zeros(3,1) A];
Bt = [D; B];
[K,~,~] = lqr(At,Bt,Q,R);
X0 = [0; x0];

[Tout,Xout]=ode45( @(t,x) myclpTrackLQ(t,x,A,B,C, D, K,0), tspan, X0 );
U=[];
for m=1:1:length(Tout)
    [xdot,u] = myclpTrackLQ(Tout(m),Xout(m,:)',A,B,C, D, K,0);
    U = [U;u];
end

% plots
figure(2);
plot(Tout,Xout(:,2),'linewidth',2);
hold on
yline(r,'linewidth',2)
legend('x_1','r')
xlabel('Time [sec]')
ylabel('x'); set(gca,'fontsize',12)
title('Problem 3D: Forward Velocity and Set-Point vs Time')

figure(3);
plot(Tout,Xout(:,1),'linewidth',2);
hold on
plot(Tout,U,'linewidth',2);
legend('x_I','u')
xlabel('Time [sec]')
ylabel('x,u'); set(gca,'fontsize',12)
title('Problem 3D: Integrator State and Control Input vs Time')
% Xout(end,2)
%% part e
w = 0.8; % disturbance additive to control input
[Tout,Xoutdist]=ode45( @(t,x) myclpTrackLQ(t,x,A,B,C, D, K,w), tspan, X0 );
Udist=[];
for m=1:1:length(Tout)
    [xdot,u] = myclpTrackLQ(Tout(m),Xoutdist(m,:)',A,B,C, D, K,w);
    Udist = [Udist;u];
end

% plots
figure(4);
plot(Tout,Xoutdist(:,2),'linewidth',2);
hold on
yline(r,'linewidth',2)
legend('x_1','r')
xlabel('Time [sec]')
ylabel('x'); set(gca,'fontsize',12)
title('Problem 3E: Forward Velocity and Set-Point vs Time with Disturbance')

figure(5);
plot(Tout,Xoutdist(:,1),'linewidth',2);
hold on
plot(Tout,Udist,'linewidth',2);
legend('x_I','u')
xlabel('Time [sec]')
ylabel('x,u'); set(gca,'fontsize',12)
title('Problem 3E: Integrator State and Control Input vs Time with Disturbance')
% Xout(end,2)
%% functions
function [xdot, u] = myclp(t,x,A,B,K)
u = K*x; % your code here
xdot = A*x+B*u; % your code here
end

function [Xtdot, u] = myclpTrackLQ(t,X,A,B,C, D, K,w)
% get the states from the augmented state vector
xIt = X(1);     % this is defined as a newly introduced state where xItdot = e
xt = X(2:end);  % my regular states
% get the feedback gains for each state from the feedback gain matrix (- of
% what lqr function outputs)
K1 = K(1); 
K2 = K(2:end);
% compute the LQ control
u = -K2*xt - K1*xIt + w;
% compute the states dot
yt = C*xt + D*u;
global r
e = yt - r;
xItdot = e; 
xtdot  = A*xt+B*u;
% output the state sifferential equations
Xtdot = [xItdot; xtdot];
end
