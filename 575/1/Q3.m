clear all; clc; close all;
% given matrices
v_min = [-55; -30; -30; -30].*pi/180;           % u = Bv
B     = [0       -4.2423   4.2423    1.4871;...
         1.6532  -1.2735  -1.2735    0.0024;...
         0       -0.2805   0.2805   -0.8823];
u_des = [0    -2.25    0.5]';

A_ineq = [ eye(4) ; -eye(4) ];
b_ineq = [-v_min  ; -v_min   ];
Q      = 2*B'*B;
% c      = -2.*u_des'*B;
c      = -2.*B'*u_des'*B;
d      = u_des'*u_des;

x = quadprog(Q,c,A_ineq,b_ineq);