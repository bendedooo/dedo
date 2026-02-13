function [statesOut] = DP(t,states,V_M,V_T)
% solve these to find the trajectory of the target
theta_T = pi/2; % bunu formulize et

beta = states(5);
beta = beta+(0.25)*randn;
R = states(6);
theta = beta; 
R_dotOut = V_T*cos(beta-theta_T) - V_M*cos(beta-theta); %ders notunde theta t yi 0 aliyo ama bizimki 90 dimi?
beta_dotOut = -(V_T*sin(beta-theta_T)-V_M*sin(beta-theta))/(R);
X_Tx_dotOut = 0;
X_Ty_dotOut = V_T;
X_Mx_dotOut = V_M*cos(beta);
X_My_dotOut = V_M*sin(beta);
statesOut = [X_Mx_dotOut; X_My_dotOut; X_Tx_dotOut; X_Ty_dotOut; beta_dotOut; R_dotOut];


% x_dot = [R_dot; beta_dot];

% x_dot = [V_T*cos(beta-theta_T) - V_M*cos(beta-theta);...
%          -(V_T*sin(beta-theta_T)-V_M*sin(beta-theta))/(R)];
end
