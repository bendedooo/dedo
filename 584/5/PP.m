function [statesOut] = PP(t,states,V_M,V_T,lambda,theta)
% solve these to find the trajectory of the target
theta_T = pi/2+ cos(t);
beta = states(5);
% betaMeas = beta+(0.25)*randn;
R = states(6);
% theta = states(7); 
R_dotOut = V_T*cos(beta-theta_T) - V_M*cos(beta-theta); %ders notunde theta t yi 0 aliyo ama bizimki 90 dimi?
beta_dotOut = -(V_T*sin(beta-theta_T)-V_M*sin(beta-theta))/(R);
X_Tx_dotOut = V_T*cos(theta_T);
X_Ty_dotOut = V_T*sin(theta_T);
X_Mx_dotOut = V_M*cos(theta);
X_My_dotOut = V_M*sin(theta);
theta_dot = lambda*beta_dotOut;
statesOut = [X_Mx_dotOut; X_My_dotOut; X_Tx_dotOut; X_Ty_dotOut; beta_dotOut; R_dotOut; theta_dot];


% x_dot = [R_dot; beta_dot];

% x_dot = [V_T*cos(beta-theta_T) - V_M*cos(beta-theta);...
%          -(V_T*sin(beta-theta_T)-V_M*sin(beta-theta))/(R)];
end
