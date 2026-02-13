function beta_dotOut = getBetaDot (states,V_T,V_M)
% take beta dot
beta = states(:,5);
R = states(:,6);
theta_T = pi/2;
theta = beta; 
beta_dotOut = -(V_T.*sin(beta-theta_T)-V_M.*sin(beta-theta))./(R);
end
