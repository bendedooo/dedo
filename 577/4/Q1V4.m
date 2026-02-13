% HW4 Q1
clear all; close all; clc
load('BurgersSolFull')
load('paramBurgerFD')

TrainRaw = solf;
timeavg = mean(TrainRaw,2);    % M = mean(A,dim) returns the mean along dimension dim. For example, if A is a matrix,
% then mean(A,2) is a column vector containing the mean of each row.
% Train = bsxfun(@minus,TrainRaw,timeavg);  % Train = Train - Time Average
% ortalamayi cikartip s

Train = TrainRaw;
[Psi,S,V] = svd(Train,'econ'); % POD modes Psi matrix

cumEnergy = cumsum(diag(S))./sum(diag(S)); % normalized cumulative energy = normalized sum of singular values
figure(2)
plot(cumEnergy, LineWidth=2);
grid on
ylabel('Cumulative Energy')
xlabel('Snapshots')
title('Cumulative Energy')

thresholdInd = find(cumEnergy>=0.9999,1);


%% solve in new coord
modeNoArray = [10:5:300];
L2errorRe   = zeros(length(modeNoArray),1);% ,size(solf,2));
L2errorPro  = zeros(length(modeNoArray),1);% ,size(solf,2));
% U_r         = zeros(size(solf,1),1);
% U_p         = zeros(size(solf,1),1);
            
u0 = solf(:,1);
dt = (tfin - t0 )/nt; tout = zeros(nt,1); t = t0; tout(1) = t0;

for ii = 1:length(modeNoArray)
    m = modeNoArray(ii);
    PhiT = Psi(:,1:m);    % Truncated POD basis
    V = PhiT;
    solr = zeros(m,nt); a = V'*u0; solr(:,1) = a;
        for jj = 1:nt-1
        t=  t + dt;
        tout(jj+1) = t;
        a = a + dt*V'*F(V*a);
        solr(:,jj+1) = a;
        end
%         % the final snapshot error
%         a_r = solr(:,end);
%         u_f = Train(:,end);
%         u_r = V*a_r;
%         L2errorRe(ii,:) = norm(u_r-u_f);
%         u_p = V*V'*u_f;
%         L2errorPro(ii,:) = norm(u_p-u_f);
        % difference of all snapshots norm
        A_r = V*solr;
        A_p = V*V'*Train;
        L2errorRe(ii,:) = norm(A_r-Train);
        L2errorPro(ii,:) = norm(A_p-Train);

%         U_r(ii,:) = u_r;
%         U_p(ii,:) = u_p;
end

figure(3)
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.1 0.7 0.7])
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_2$', 'Interpreter','latex' )
title('ROM Reconstruction Error vs Dimension')

figure(4)
semilogy(modeNoArray,L2errorPro,'LineWidth',2,Color=[0.1 0.1 0.7])
hold on
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.1 0.7 0.7])
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_2$', 'Interpreter','latex' )
title('ROM Reconstruction and Projection Error vs Dimension')
legend('Projection Error','Reconstruction Error')

save('hw4','L2errorRe','modeNoArray')
%% part c

load('BurgersSolFullc')
load('paramBurgerFDc') 
dt = (tfin - t0 )/nt; tout = zeros(nt,1); t = t0; tout(1) = t0;
Train = solf;
l = ceil(length(Train)/2);
for ii = 1:length(modeNoArray)
    m = modeNoArray(ii);
    PhiT = Psi(:,1:m);    % Truncated POD basis
    V = PhiT;
    solr = zeros(m,nt); a = V'*u0; solr(:,1) = a;
        for jj = 1:nt-1
        t=  t + dt;
        tout(jj+1) = t;
        aaa = F(V*a);
        a = a + dt*V'*F(V*a);
        solr(:,jj+1) = a;
        end
%         a_r = solr(:,end);
%         u_f = Train(:,end);
%         u_r = V*a_r;
%         L2errorRe(ii,:) = norm(u_r-u_f);
%         u_p = V*V'*u_f;
%         L2errorPro(ii,:) = norm(u_p-u_f);
%         U_r(ii,:) = u_r;
%         U_p(ii,:) = u_p;
        A_r = V*solr;
        A_p = V*V'*Train;
        L2errorRe(ii,:) = norm(A_r(l:end)-Train(l:end),'fro');
%         L2errorPro(ii,:) = norm(A_p-Train);
end
figure(5)
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.1 0.7 0.7])
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_2$', 'Interpreter','latex' )
title('ROM Reconstruction Error vs Dimension')

figure(6)
semilogy(modeNoArray,L2errorPro,'LineWidth',2,Color=[0.1 0.1 0.7])
hold on
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.1 0.7 0.7])
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_2$', 'Interpreter','latex' )
title('ROM Reconstruction and Projection Error vs Dimension')
legend('Projection Error','Reconstruction Error')
%%
solf1  = [zeros(nt,1) A_r' zeros(nt,1)]; % add boundaries
xx    = linspace(x0,xf,n+2);

figure(1)
surfc(xx,tspan ,solf1);
shading interp
title(['Sol of Full System (FD):dim = ' num2str(n)]);
xlabel('x');
ylabel('t');
zlabel('u');