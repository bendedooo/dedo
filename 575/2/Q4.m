clear all; clc; close all;
 
W = 75000;         %[lbs]
S = 950;           %[ft^2]
AR = 5.9;
C_D_0 = 0.015;
K = 0.05;
rho = 0.00238;     %[slug/ft^3]
C_L_alpha = 0.12;
T_max = 11000;     %[lbs]
C_L_max = 2.8;

syms V q1 q2
 % function to be minimized
f  = V;           
% inequality constraints <= 0
h1 = 2*W/(rho*S*V^2) - C_L_max;
h2 = 0.5*rho*V^2*S*C_D_0 + 2*K*W^2/(rho*V^2*S) - T_max; 

df = diff(f,V);
dh1 = diff(h1,V);
dh2 = diff(h2,V);

q0 = 1;                             % KKT conditions are assumed to hold, q0=1
dL = q0*df + q1*dh1  + q2*dh2;      % gradient of the lagrangian

% h1 is active
[V1cand] = solve(h1 == 0, V);   % solve for V in the active constraint
V1cand   = double(V1cand(V1cand>0))
V1 = V1cand;
h2V1 = double(subs(dh1, [V], [V1]));
V1 = sqrt(2*W/(rho*S*C_L_max)); %bu neden 0 vermiyo


% h2 is active
[V2cand] = solve(h2 == 0, V);   % solve for V in the active constraint
V2cand   = double(V2cand(V2cand>0))


% CASES 
A={1} % solve for q1
% dh1d = double(subs(dh1, [V q2], [V1 0]));    % obtain dh1 value with that V value
dL1 = (subs(dL, [V q2], [ V1 0]));    % obtain dh1 value with that V value
V1
[q1Sol] = double(solve(dL1 == 0, q1))              % solve for q1 in the gradient of the lagrangian

A={2} % solve for q2
for i = 1:length(V2cand)   %check for which V values I get q>0
    V2   = V2cand(i);
h1V2 = double(subs(dh1, [V], [V2]));

dL2 = (subs(dL, [V q1], [ V2 0]));    % obtain dh1 value with that V value
[q2Sol] = double(solve(dL2 == 0, q2));              % solve for q1 in the gradient of the lagrangian
if q2Sol>0
    V2
    q2Sol
end
end
A={1,2}
