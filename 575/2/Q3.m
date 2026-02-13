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

syms V q1

f  = -V;           % function to be minimized

% inequality constraints
h1 = 0.5*rho*V^2*S*C_D_0 + 2*K*W^2/(rho*V^2*S) - T_max; % <= 0

% L = q0*f + q1*h1;
% dL = q0*df + q1*dh1;

df = diff(f,V);
dh1 = diff(h1,V);

q0 = 1;                                 % KKT conditions are assumed to hold, q0=1

[V1cand] = solve(h1 == 0, V);           % solve for V in the active constraint
V1cand   = double(V1cand(V1cand>0));
V1   = max(V1cand);

for i = 1:length(V1cand)                % check for which V values I get q>0 (dual feasibility constraint )
    V1   = V1cand(i);
    dh1 = double(subs(dh1, [V], [V1]));     % obtain dh1 value with that V value
    dL = q0*df + q1*dh1;                    % gradient of the lagrangian
    [q1Sol] = double(solve(dL == 0, q1));   % solve for q1 in the gradient of the lagrangian (stationarity constraint)
    if q1Sol>0                              % dual feasibility constraint
        V1
        q1Sol
    end
end
