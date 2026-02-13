%
%    n =  inner grid points (not include the boundary points)
%
% OUTPUT
%
%   solf   = n-by-nt matrix of the solution at tj, j = 1,..., nt
%   t_full = CPU time for solving full-order system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Nov 4, 2010
%

n  = 800
x0 = 0; 
xf = 4;
t0 = 0;
tfin = 2.5;%1.5*2;

%coeff matrices from FD  discretization

dx = (xf- x0)/(n+1)
d1 = repmat(1/dx,n,1);
d2 = repmat(1/(dx^2),n,1);
Ax = spdiags([-d1 d1  0*d1], -1:1, n, n); 
F = @(u) -0.5*Ax*(u.*u);

% param for solving ODEs

xs  = linspace(x0, xf, n+2)';
u0 = initial_cond(xs(2:n+1));

dt = 0.2*dx
nt     = floor(tfin/dt); 
tspan  = linspace(t0,tfin,nt);

save paramBurgerFDc.mat x0 xf t0 tfin n dx Ax F u0 nt tspan

%% Solve using FD 
%%%%%%% solve by semi-implicite Euler %%%%%%%%%%%%%%%%%%%

dt = (tfin - t0 )/nt

tout = zeros(nt,1);
tout(1) = t0; 
t = t0;

solf = zeros(n,nt);
solf(:,1) = u0; 
u = u0;

tic
for j = 1:nt-1,
    
    t=  t + dt; 
    tout(j+1) = t;  
    u = u + dt*F(u) ;  
    solf(:,j+1) = u;
    
end
t_full = toc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save BurgersSolFullc.mat tout solf


%% Plot
solf1  = [zeros(nt,1) solf' zeros(nt,1)]; % add boundaries
xx    = linspace(x0,xf,n+2);

figure(1)
surfc(xx,tspan ,solf1);
shading interp
title(['Sol of Full System (FD):dim = ' num2str(n)]);
xlabel('x');
ylabel('t');
zlabel('u');
