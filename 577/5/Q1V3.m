% HW4 Q1
clear all; close all; clc
load('BurgersSolFull')
load('paramBurgerFD')
load('hw4.mat');

% Compute the POD modes
Train     = solf(:,1:end-1);

[Psi,S,V] = svd(Train,'econ');

dt = (tfin - t0 )/nt;
% solNL = F(Train);
solNL = Train(:,2:end);

[PsiNL,~,~] = svd(solNL,'econ');
% PsiNL = PhiT;


for ii = 1:length(modeNoArray)
    m = modeNoArray(ii);
    pointNo = m; % for QDEIM

    PhiT = Psi(:,1:m);    % Truncated POD basis
    PsiNLT = PsiNL(:,1:m);

%     PsiNLT = PhiT;
    [~,~,pivot] = qr(PsiNLT','vector');
    sensors = pivot(1:pointNo);
    
    P = zeros(m,800);
    for b=1:m
        P(b,sensors(b))=1;
    end


    % compute the offline matrix
    offlineMat = PhiT' * PsiNLT * pinv(PsiNLT(sensors,:));

    u0 = Train(:,1);
    %     u0 = solNL(:,1);
    a = PhiT'*u0;
    solr = zeros(m,nt);
    solr(:,1) = PhiT'*u0; % set the initial conditions

    for jj = 1:nt-1
        a = a + dt*offlineMat*P*F(PhiT*a);
        solr(:,jj+1) = a;
    end

    %     % the final snapshot error
    %     a_r = solr(:,end);
    %     u_f = Train(:,end);
    %     u_r = PhiT*a_r;
    %     L2errorRe(ii,:) = norm(u_r-u_f);

    % difference of all snapshots norm
    A_r = PhiT*solr;
    L2errorReM(ii,:) = norm(A_r-solf,'fro');
end

figure(1)
semilogy(modeNoArray,L2errorReM,'LineWidth',2,Color=[0.1 0.7 0.7])
hold on
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.3 0.3 0.7])
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_2$', 'Interpreter','latex' )
title('ROM Reconstruction Error vs Dimension')
legend('with QDEIM','without QDEIM',Location='best')



% %%
% solf1  = [zeros(nt,1) A_r' zeros(nt,1)]; % add boundaries
% xx    = linspace(x0,xf,n+2);
%
% figure(1)
% surfc(xx,tspan ,solf1);
% shading interp
% title(['Sol of Full System (FD):dim = ' num2str(n)]);
% xlabel('x');
% ylabel('t');
% zlabel('u');