function [X22, P22, C_A2B2] = KFP3AV2(X11,C_A2B1,P11,accMeas,omegaMeas,posMeas,measAvailable,T)
g = 9.80665;
gVec = [0 0 -g]'; % in frame A

% KF Matrices
Q = 10*eye(6);
R = 0.001*eye(3); % meas noise
A = [eye(3) T*eye(3);
    zeros(3) eye(3)];
B = [T^2/2*eye(3);
     T*eye(3)];

% Poissons Eqn
% Xi = expm(-T*omegaSS(omegaMeas));
nk = omegaMeas/norm(omegaMeas);
Xi = expm(-norm(omegaMeas)*T*omegaSS(nk)); % OMEGAk
% Xi = -omegaSS(omegaMeas); % OMEGAk


% Double Integration
u = C_A2B1'*accMeas - gVec; % ?????? C K+1
X21 = A*X11+B*u;

% K
C = [eye(3) zeros(3)]*measAvailable;
P21 = A*P11*A'+Q;
K = P21*C'*inv(C*P21*C'+R);

% posterior
X22 = X21 + K*(posMeas-C*X21);
C_A2B2 = Xi*C_A2B1;
P22 = P21-K*C*P21;

end