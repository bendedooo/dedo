% HW 5 P6: Noisy DP and PP
clc
clear all
close all


p1_0 = 1;
p2_0 = 0;
timeVec = flip(0:0.001:9.999);
states0 = [p1_0; p2_0];
states0 = [p1_0; p2_0];

lambdaArray = [3 4 5];
lambda = 3;
for i=1:3
    lambda = lambdaArray(i);

    odeFunc = @(t,x)Adj(t,x,lambda);
    optPos = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    [tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);

    nTG = -lambda*statesOut(:,2)./(10-tOut);
    figure(1)
    plot(10-tOut,nTG,'LineWidth',2)
    hold on
    

end

xlabel('$\frac{\bar {t_f},t)}{T}\  time-to-go [s]$','interpreter','latex');
ylabel('$-TG(\bar {t_f},t) [m]$','interpreter','latex')
grid on
legend('\Lambda = 3','\Lambda = 4','\Lambda = 5',...
    'Location','best')
title('Problem 7')
