% HW4 Q1
clear all; close all; clc
load('BurgersSolFull')
load('paramBurgerFD')
load('hw4.mat');

% Compute the POD modes
Train     = solf;
[Psi,S,V] = svd(Train,'econ');

dt = (tfin - t0 )/nt;
% solNL = F(Train);
% solNL = Train(:,2:end);
% [PsiNL,~,~] = svd(solNL,'econ');
% PsiNL = PhiT;
% [~,~,pivot] = qr(PsiNL','vector');


for ii = 1:length(modeNoArray)
    m = modeNoArray(ii);
    PhiT = Psi(:,1:m);    % Truncated POD basis
%     PsiNLT = PsiNL(:,1:m);
    PsiNLT = PhiT;
    [~,~,pivot] = qr(PsiNLT','vector');

    pointNo = m; % for QDEIM
    sensors = pivot(1:pointNo);

    P = zeros(m,800);
    for jj=1:m
        P(jj,sensors(jj))=1;
    end

    % compute the offline matrix
    offlineMat = PhiT' * PsiNLT * pinv(PsiNLT(sensors,:));

    u0 = Train(:,1);
    %     u0 = solNL(:,1);
    a = PhiT'*u0;
    solr = zeros(m,nt);
    solr(:,1) = PhiT'*u0; % set the initial conditions

    for jj = 1:nt-1
        a = a + dt*offlineMat*P*Fu(PhiT*a);
        solr(:,jj+1) = a;
    end

    % difference of all snapshots norm
    A_r = PhiT*solr;
    L2errorReM(ii,:) = norm(A_r-Train,'fro');
end

figure(1)
semilogy(modeNoArray,L2errorReM,'LineWidth',2,Color=[0.1 0.7 0.7])
hold on
semilogy(modeNoArray,L2errorRe,'LineWidth',2,Color=[0.3 0.3 0.7])
grid on
xlabel('POD-Galerkin ROM Dimension')
ylabel('$|\!|\ \tilde{u}-u|\!|_{Fro}$', 'Interpreter','latex' )
title('ROM Reconstruction Error vs Dimension')
legend('with QDEIM','without QDEIM',Location='best')
