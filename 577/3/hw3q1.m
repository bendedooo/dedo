a1 = [5;2;6;2];
v1 = a1;
u1 = a1./norm(a1);
a2 = [6;2;6;3];
v2 = a2-dot(a2,u1)*u1;
u2 = v2/norm(v2);
a3 = [6;2;2;6];
v3 = a3-dot(a3,u2)*u2-dot(a3,u1)*u1;
u3 = v3/norm(v3);
a4 = [8;8;8;7];
v4 = a4-dot(a4,u3)*u3-dot(a4,u2)*u2-dot(a4,u1)*u1;
u4 = v4/norm(v4);
A  = [a1 a2 a3 a4];
Q  = [u1 u2 u3 u4]
R  = Q'*A
[QM,RM] = qr(A)
[Q,R,P] = qr(A)


% [QM,RM] = qr(A);
% 
% Q = [];
% A = [a4 a1 a2 a3];
% v1 = a4;
% u1 = a4./norm(a4)
% Q1 = [u1];
% R1 = Q1'*A;
% A = [a4 a2 a1 a3];
% v2 = a2-dot(u1,a2)*u1;
% u2 = v2/norm(v2);
% Q2 = [u2];
% R2 = Q2'*A;
% A = [a4 a2 a3 a1];
% v3 = a3-dot(u2,a3)*u2-dot(u1,a3)*u1;
% u3 = v3/norm(v3);
% Q3 = [u3];
% R3 = Q3'*A;
% v4 = a1-dot(u3,a1)*u3-dot(u2,a1)*u2-dot(u1,a1)*u1;
% u4 = v4/norm(v4);
% Q4 = [u4];
% R4 = Q4'*A;
% 
% Q = [Q1 Q2 Q3 Q4]
% R = [R1; R2; R3; R4]

