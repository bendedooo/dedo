function [X22, P22] = KFP2(state11,P11,L,Y)
R2 = 0.1*eye(3);
A = eye(size(state11,1));
state21 = A*state11;
P21 = P11;
C2 = CMatrix(state11,L);
K2 = P21*C2'*inv(C2*P21*C2'+R2);
P22 = P21-K2*C2*P21;
g = measModel(state21,L);
X22 = state21 + K2*(Y-g);
end