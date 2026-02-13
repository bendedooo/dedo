%% Q2
clc; clear all; close all;
load('iss'); % vector of the Hankel singular values hsv, a frequency vector w and the corresponding frequency response mag.
% create the system impulse response
A = full(A);
B = full(B);
C = full(C);
C = C(1,:);
D = 0;
sys = ss(A, B, C, D);
evalFOM = eig(A,'vector');

r = 20;  % number of retained modes
%Compute the direct system impulse response given the system and the time vector.
tf = 0.1;
dt = 5e-05;
timeVec = 0:dt:tf;

sysD = c2d(sys,dt);
evalFOMD = eig(sysD.A,'vector');



[yFullA,t,xFullA] = impulse(sys,timeVec); % impulse assumes initial state values are zero
% xFull: size m (# of snapshots) x n (system dim) x p/q (# of inputs/outputs)

%???---POD boyle mi yapiyoruz yoksa---???
%% Part a: POD Galerkin
%  
% olhpde degil unstable 
X1=[]; X2=[];
for i=2:size(xFullA,1) % Start at 2 to avoid the D matrix
%       XSnapshots = [XSnapshots squeeze(xFullA(i,:,:))];
      X1 = [X1 squeeze(xFullA(i,:,:))];
%       X2 = [X2 squeeze(xFullA(i+1,:,:))];
end
[Psi,S,V] = svd(X1,'econ'); % POD modes Psi matrix
PsiT = Psi(:,1:r);
%?????

% AGI = X2*pinv(X1);
% sysD = c2d(sys,dt);
AGI  = sys.A;

F = @(u) PsiT*PsiT'*AGI*(u);
AGal = PsiT'*AGI*PsiT;

solr =[];

    u0 = squeeze(xFullA(1,:,1:3));
    u0 = X1(:,1:3);
    
    a = PsiT'*u0; solr = a;
    for jj = 1:size(xFullA,1)-1
        a = a + dt*PsiT'*F(PsiT*a);
        solr = [solr a];
    end
  
    % the final snapshot error
%     a_r = solr(:,end-3:end);
%     u_f = X1(:,end);
%     u_r = PsiT*a_r;
%     L2errorRe(ii,:) = norm(u_r-u_f);
    % difference of all snapshots norm
    A_r = PsiT*solr;
%     L2errorReM(ii,:) = norm(A_r-X1);
    % matrix norm sikinti cikariyo

% X1G=A_r(:,1:end-3); 
% X2G=A_r(:,4:end);
% AGO = X2G*pinv(X1G);


X1G=solr(:,1:end-3); 
X2G=solr(:,4:end);
AGO = X2G*pinv(X1G);

    evalGal = eig(AGal,'vector');
    evalFOM = eig(A,'vector');

    figNo = 1;
    figure(figNo)
    plot(evalGal,'b*');
    hold on
    plot(evalFOM,'r*');
    grid on
    ylabel('Im')
    xlabel('Re')
    title('Eigenvalue Spectrum ERA')
    legend('ROM','FOM')

    reEvalGal = real(evalGal)>=0; % eigval will be all zeros if eigvals are on the OLHP
    any(reEvalGal) % B = any(A) returns logical 1 ( true ) if any of the elements of A is a nonzero number or is logical 1
        % and returns logical 0 ( false ) if all the elements are zero
        % if this value is 0 -> stable all eigvals on OLHP
        % if this value is 1 -> stable all eigvals on OLHP

%% Part b: ERA
% bu farkli sadece simdilik
% discrete system check the unit circle
tfArray = [0.1];%, 0.2, 0.5
H =[];
H2 =[];

% ------ bu tf artircaz mi azaltcaz mi??? ------- %
for cc = 1:length(tfArray)
    tf      = tfArray(cc);
    timeVec = 0:dt:tf;
    Y = permute(yFullA,[2 3 1]);
    [q, p, m] = size(Y); % Y: q x p x m; q: number of outputs
                                       % p: number of inputs
                                       % m: number of snapshots

    % Hankel Matrix q(mo) x p(mp); mp and mo are hyperparameters, mo=mp=m/2
    mo = m/2;
    mp = m/2;
%     Dr = []; 
    M = [];

    for ii = 1:q
        for jj = 1:p
%             Dr(ii,jj) = Y(ii,jj,1);        % 1st block of Y: D matrix
            M(ii,jj,:) = Y(ii,jj,2:end);   % romove D to build the Hankel matrix
        end
    end

    for ii = 1:mo
        for jj = 1:mp
            for kk = 1:q
                for ll = 1:p
                    H (q*ii-q+kk, p*jj-p+ll) = M(kk, ll, ii+jj-1);  % Hankel
                    H2(q*ii-q+kk, p*jj-p+ll) = M(kk, ll, ii+jj);    % Shifted Hankel
                end
            end
        end
    end

    [U,S,V] = svd(H,'econ');
    Sigma   = S(1:r,1:r);
    Ur      = U(:,1:r);
    Vr      = V(:,1:r);
    Ar      = Sigma^(-0.5)*Ur'*H2*Vr*Sigma^(-0.5);
    Br      = Sigma^(-0.5)*Ur'*H2(:,1:p);
    Cr      = H(1:q,:)*Vr*Sigma^(-0.5);
    Dr      = 0;

    sysERA = ss(Ar, Br, Cr, Dr);
    evalERA = eig(Ar,'vector');


    figNo = cc+1;
    figure(figNo)
    plot(evalERA,'b*');
    hold on
    plot(evalFOMD,'r*');
    grid on
    ylabel('Im')
    xlabel('Re')
    title('Eigenvalue Spectrum ERA')
    legend('ROM','FOM')

    reEvalERA = real(evalERA)>=0; % eigval will be all zeros if eigvals are on the OLHP
    any(reEvalERA) % B = any(A) returns logical 1 ( true ) if any of the elements of A is a nonzero number or is logical 1
        % and returns logical 0 ( false ) if all the elements are zero
        % if this value is 0 -> stable all eigvals on OLHP
        % if this value is 1 -> stable all eigvals on OLHP
end
%% Part c: DMD
% bu ayni
% Y1 = squeeze(permute(yFull(1:end-1,:,:),[2 3 1]));
% Y2 = squeeze(permute(yFull(2:end  ,:,:),[2 3 1]));
% 
% %---- bu y1 ve y2 snapshot matrixi hankel gibi olusturmucaz dimi??? ----%
% %---- bu svd 3 tane cikiyo ----%
%     X = permute(xFull,[2 3 1]);
%     [q, p, m] = size(X); % X: q x p x m; q: number of states
%                                        % p: number of inputs
%                                        % m: number of snapshots
%     X1 = [squeeze(X(:,1,1:end-1)); squeeze(X(:,2,1:end-1)); squeeze(X(:,3,1:end-1))];



% Burda Xleri mi kullancaz???

X1 = []; X2 = [];
for i=2:size(xFullA,1)-1 % Start at 2 to avoid the D matrix
      X1 = [X1 squeeze(xFullA(i,:,:))];
      X2 = [X2 squeeze(xFullA(i+1,:,:))];
end

% r = rank(X1);

[U,S,V] = svd(X1,'econ');
Sigma   = S(1:r,1:r);
Ur      = U(:,1:r);
Vr      = V(:,1:r);
Ar      = Ur'*X2*Vr*inv(Sigma);
[W,evalDMD]   = eig(Ar,'vector'); % AW=WD
Phi     = X2*Vr*inv(Sigma)*W;
% Ard = c2d(Ar);
% eigenvaluelari ayni olmali
figure(3)
plot(evalDMD,'b*');
hold on
plot(evalFOMD,'r*');
grid on
ylabel('Im')
xlabel('Re')
title('Eigenvaluee Spectrum DMD')
legend('ROM','FOM')

reEvalDMD = real(evalDMD)>=0; % eigval will be all zeros if eigvals are on the OLHP
any(reEvalDMD) % B = any(A) returns logical 1 ( true ) if any of the elements of A is a nonzero number or is logical 1
% and returns logical 0 ( false ) if all the elements are zero
% if this value is 0 -> stable all eigvals on OLHP
% if this value is 1 -> stable all eigvals on OLHP

%% Part d: POD_Galerkin vs ERA vs BT
rArray = [5, 10, 20, 40];

for ii = 1:length(rArray)
    r = rArray(ii)

    %% POD_Galerkin




    ePGAL = eig(Ar(1:r,1:r),'vector');
    reEvalPGAL = real(ePGAL)>=0; % eigval will be all zeros if eigvals are on the OLHP
    isUnstablePGAL = any(reEvalPGAL)


    %% ERA

    eERA = eig(Ar(1:r,1:r),'vector');
    reEvalERA = real(eERA)>=0; % eigval will be all zeros if eigvals are on the OLHP
    isUnstableERA = any(reEvalERA)

    %% BT
    % ----- sistemi biliyo muyuz yani ?? -------
    Wc = gram(sys, 'c');
    Wo = gram(sys, 'o');
    %Compute the unscaledtransformation matrix.
    [Tu, D] = eig(Wc*Wo);
    %Transform the system to the new coordinates.
    ATu = inv(Tu) * A * Tu;
    BTu = inv(Tu) * B;
    CTu = C * Tu;
    DTu = 0;
    syst = ss(ATu, BTu, CTu, DTu);
    %Compute the Gramians of the transformed system.
    Sigma_c = gram(syst, 'c');  %These gramians are diagonal but notequal
    Sigma_o = gram(syst, 'o');
    %Compute the scaling factor.
    Sigma_s = diag(Sigma_c) ./ diag(Sigma_o);
    %Compute the balancing transformation.
    T = Tu * diag(Sigma_s .^ (1/4));
    %Compute the balanced Gramians (diagonals of Sigma are the Hankel singular values).
    Sigma = diag(Sigma_c) .^ (1/2) .* diag(Sigma_o) .^ (1/2);
    %Sort the singular values (and their correponding modes) in descending order.
    %In larger systems, here you can truncate lower modes to generate a "reduced-order" model.
    [Sigma, permind] = sort(Sigma, 'descend');
    T = T(:, permind);
    %Compute the balanced ROM matrices.
    A_b = inv(T) * A * T;
    B_b = inv(T) * B;
    C_b = C * T;
    D_b = 0;
    sysb = ss(A_b, B_b, C_b, D_b);
    %test
    Wc_new = gram(sysb, 'c');
    Wo_new = gram(sysb, 'o');
    eBT = eig(A_b(1:r,1:r),'vector');
    reEvalBT = real(eBT)>=0; % eigval will be all zeros if eigvals are on the OLHP
    isUnstableBT = any(reEvalBT)








    
    
    
    
    
end
    
    %% cop








    HankelOC = []; % Compute Hankel matrix H=OC
    for i=2:size(xAdj,1) % Start at 2 to avoid the D matrix
        Hrow = [];
        for j=2:size(xFull,1)
            Ystar = permute(squeeze(xAdj(i,:,:)),[2 1]);
            MarkovParameter = Ystar*squeeze(xFull(j,:,:));
            Hrow = [Hrow MarkovParameter];
        end
        HankelOC = [HankelOC; Hrow];
    end
    [U,Sig,V] = svd(HankelOC);




    %???---POD boyle mi yapiyoruz yoksa---???
    % C = Train'*Train; [V,D] = eig(C); Sigma =
    %diag(sqrt(D)); Phi = Train*V*inv(Sigma);