function [X2, C_A2B2] = KFP3(X1,C_A2B1,Y,T)
g = 9.80665;
gVec = [0 0 -g]'; % in frame A
Y1 = Y;
omegaMeas = Y(1:3);
accMeas = Y(4:6);
nk = omegaMeas/norm(omegaMeas);
OMEGAk = exp(-norm(omegaMeas)*T*omegaSS(nk));
Y2 = OMEGAk*Y1;

% PE
Ak = expm(-norm(omegaMeas)*T*(nk));
att21 = Ak*att11;
att22 = att21; + K*0;

% DI
pos21 = Ak*posn11+B*gVec+Bk*Y1;
posn22 = posn21 + K*0;

Ad = [eye(3) T*eye(3);
     zeros(3) eye(3)];
Bd = [T^2/2*eye(3);
     T*eye(3)];
u = Y1'*accMeas - gVec; % ?????? C K+1
X2 = Ad*X1 + Bd*u;

%% Nav
A = [eye(3) T*eye(3);
     zeros(3) eye(3)];
B = [T^2/2*eye(3);
     T*eye(3)];
% C = ;
% D = ;

omegaMeas = Y(1:3);
C_A2B2 = (exp(-T.*(-omegaSS(omegaMeas))))*C_A2B1;
accMeas = Y(4:6);
u = C_A2B1'*accMeas - gVec; % ?????? C K+1
X2 = A*X1 + B*u;
% 
% Q = ;
% R = 0.1*eye(3); % meas noise
% 
% P21 = A*P11*A'+Q;
% K = P21*C2'*inv(C2*P21*C2'+R);
% X21 = A*X11 + B*u;
% 
% % posteriors
% P22 = P21-K*C2*P21;
% X22 = X21 + K(Y2-C2*X21);



end