function beta_ddotOut = getBetaDDot (states,V_T,V_M,beta_dot)
% take beta dot
beta = states(:,5);
R = states(:,6);
theta_T = pi/2;
theta = beta; 
R_dot = V_T.*cos(beta-theta_T) - V_M.*cos(beta-theta); 
beta_ddotOut = (-R.*V_T.*cos(beta-theta_T).*beta_dot-(V_M.*sin(beta-theta)-V_T.*sin(beta-theta_T).*R_dot))./(R.^2);
end