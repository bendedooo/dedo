% Q4 Part c
clear all; close all; clc;
J = 1;
B = 1/J;
% x: theta dot
N = 20;
dT = 0.25;
% terminal constraints
theta_f = 0;
theta_d_f = 0;
xN = [theta_f theta_d_f]';

x0 = [pi/2, 0]';

R = 1;
Q = 0;

A = [0 1; 0 0];
B = [0; 1/J];


% AA = [A B;...
%  zeros(1,length(A)+1) ]
% XX = expm(AA*dT)
% Ad = 
Ad = expm(A*dT)
Bd = (dT*eye(length(A))+A.*(dT^2/2))*B


PHI = [ Ad+Bd*inv(R)*Bd'*inv(Ad)'*Q ,    -Bd*inv(R)*Bd'*inv(Ad)';...
        -Ad'*Q                 ,    inv(Ad)'                 ];
G   = PHI^N;
% p0 = inv(G(1,1))*(xN - G(1,1)*x0);

% [x0; p0]= inv(PHI)*[xN; pN]

p = inv([G(1:2,3:4) zeros(2) ;G(3:4,3:4) -eye(2)])*[-G(1:2,1:2)*x0+xN;-G(3:4,1:2)*x0];
p0 = p(1:2);
pN = p(3:4);
% pk = p0;
% xk = x0;

X = [x0];
P = [p0];
U = [];

for ii=1:1:N-1
    P(:,ii+1) = inv(Ad)'*(P(:,ii) - Q*X(:,ii));
    U(:,ii) = -inv(R)*Bd'*P(:,ii+1);
    X(:,ii+1) = Ad*X(:,ii) + Bd*U(:,ii);
    
end

% P(ii+1) = pN;
X(:,ii+1) = xN;
%%
figure(1)
plot(X(1,:),X(2,:),'LineWidth',2)
hold on
scatter(X(1,1),X(2,1),'filled','r')
xlabel('x_1'); ylabel('x_2'); grid on;
title('Problem 4C Trajectory')
legend('Satellite Trajectory','Spacecraft Initial Point')

figure(2)
subplot(2,1,1)
plot(P(1,:),'LineWidth',2)
xlabel('k'); ylabel('p_1'); grid on;
title('Problem 4C Adjoint Variable 1 vs step')
subplot(2,1,2)
plot(P(2,:),'LineWidth',2)
xlabel('k'); ylabel('p_2'); grid on;
title('Problem 4C Adjoint Variable 2 vs step')

figure(3)
plot(U,'LineWidth',2)
xlabel('k'); ylabel('u'); grid on;
title('Problem 4C Control Input vs Time')
