function out=gFunc2(t,a,dd,V,Phi,dx)
a0 = [1 a']';               % initial condition for Galerkin POD ODE 

gIn = Phi*a0;
Phidot = [zeros(1,128); (Phi(1:end-1,:)-Phi(2:end,:))/dx]; % boyutu tutmuyorrrrr
g_u=  -(Phi.*Phidot)*a0.^2; % a^t phi phidot ones(800,1) [1x800] mu olmalii????
out = V'*g_u;

-0.5*Ax*(u.*u);


end