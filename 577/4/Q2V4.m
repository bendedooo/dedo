%% Q2
clc; clear all; close all;
load('iss'); % vector of the Hankel singular values hsv, a frequency vector w and the corresponding frequency response mag.
% Part 1
A = full(A);
B = full(B);
C = full(C);
D = 0;

%???? dogru datayi mi kullanyorum????

sys= ss(A, B, C, D);
Wc = gram(sys, 'c'); % or Gc = S'*S;
Wo = gram(sys, 'o'); % or Go = R'*R;

%Compute the unscaledtransformation matrix.
[Tu, D] = eig(Wc*Wo);

%Transform the system to the new coordinates.
ATu = inv(Tu) * A * Tu;
BTu = inv(Tu) * B;
CTu = C * Tu;
DTu = 0;
syst = ss(ATu, BTu, CTu, DTu);

%Compute the Gramians of the transformed system.
Sigma_c = gram(syst, 'c');  %These gramians are diagonal but not equal.
Sigma_o = gram(syst, 'o');

%Compute the scaling factor.
Sigma_s = diag(Sigma_c) ./ diag(Sigma_o);

%Compute the balancing transformation.
T = Tu * diag(Sigma_s .^ (1/4)); % untruncated T

%Compute the balanced Gramians (diagonals of Sigma are the Hankel singular values).
Sigma = diag(Sigma_c) .^ (1/2) .* diag(Sigma_o) .^ (1/2);
hxvOut = Sigma;
[aa, bb] = eig(Sigma_c*Sigma_o); % U
[aaa, bbb] = eig(adjoint(Sigma_c*Sigma_o)); % V
HHH = aa * sqrt(bb) * aaa'; % H = U*Sigma*V'
cc = sqrt(diag(bb));
[UH,SH,VH] = svd(HHH,"vector");
SING = sort(SH,'descend'); % so they match
ff = Sigma-SING;
% plot(ff)

%Sort the singular values (and their correponding modes) in descending order.
%In larger systems, here you can truncate lower modes to generate a "reduced-order" model.
[Sigma, permind] = sort(Sigma, 'descend');

figure(1)
plot(Sigma,'.');%,'color',.7*[1 1 1]);
set(gca,'yscale','log','xlim',[0 length(Sigma)+1]);
grid on
% set(gca,'ylim',[1 10^5]);
ylabel('\sigma_{snapshot}')
xlabel('Snapshots')
title('Singular Value Spectrum')


cumEnergy = cumsum(Sigma)./sum(Sigma); % normalized cumulative energy = normalized sum of singular values
figure(2)
plot(cumEnergy, LineWidth=2);
set(gca,'yscale','log','xlim',[0 length(Sigma)+1]);
grid on
ylabel('Cumulative Energy')
xlabel('Snapshots')
title('Cumulative Energy')

thresholdInd = find(cumEnergy>=0.9999,1);
%% Part 2
thresholdIndArray = [thresholdInd, 5, 20, 30, 40];
for ii = 1:length(thresholdIndArray)
    thresholdInd = thresholdIndArray(ii)
    %Compute the balanced ROM matrices.
    A_b = pinv(T) * A * T;
    B_b = pinv(T) * B;
    C_b = C * T; 
    D_b = 0;
    sysb = ss(A_b, B_b, C_b, D_b);
    % %test
    % Wc_new = gram(sysb, 'c');
    % Wo_new = gram(sysb, 'o');

    eROM = eig(A_b(1:thresholdInd,1:thresholdInd),'vector');
    eFOM = eig(full(A),'vector');

    figure(ii+2)
    plot(eROM,'b*');
    hold on
    plot(eFOM,'r*');
    % set(gca,'yscale','log','xlim',[0 length(Sigma)+1]);
    grid on
    % set(gca,'ylim',[1 10^5]);
    ylabel('Im')
    xlabel('Re')
    title('Eigenvalues')
    legend('ROM','FOM')

    eigval = real(eROM)>=0; % eigval will be all zeros if eigvals are on the OLHP
    any(eigval) % B = any(A) returns logical 1 ( true ) if any of the elements of A is a nonzero number or is logical 1
    % and returns logical 0 ( false ) if all the elements are zero
    % if this value is 0 -> stable all eigvals on OLHP
    % if this value is 1 -> stable all eigvals on OLHP

end