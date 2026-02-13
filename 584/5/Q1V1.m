% HW 5 P1: Direct Pursuit
clc
clear all
close all

X_M0  = [0 0]';  % pusuer, XM notes
X_M0x = X_M0(1);
X_M0y = X_M0(2);

X_T0  = [100 0]';  % evader, XT notes
X_T0x = X_T0(1);
X_T0y = X_T0(2);

V_T = 5;
% theta_T = pi/2; % bunu formulize et

V_MArray= [6 8 11]; % for 8 (gamma>3/2), for 11 (gamma>2)

R0 = sqrt((X_T0x-X_M0x)^2+(X_T0y-X_M0y)^2); % Range
beta0 = atan2((X_T0y-X_M0y),(X_T0x-X_M0x));
states0 = [X_M0; X_T0; beta0; R0];
timeVec = 0:0.01:100;

statesOutMatrix = [];

for i = 1:length(V_MArray)
    V_M = V_MArray(i);

    odeFunc = @(t,x)Rdotbetadot(t,x,V_M,V_T);


    optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
    [tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);


%     X_TOutx = statesOut(:,3);
%     X_TOuty = statesOut(:,4);
%     X_MOutx = statesOut(:,1);
%     X_MOuty = statesOut(:,2);


%     figure(i)
%     scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
%     hold on
%     scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
%     plot(X_TOutx, X_TOuty,'LineWidth',2,'Color',[.5 0 .5])
%     plot(X_MOutx', X_MOuty','LineWidth',2,'Color',[0, 0.5, 0.5])
%     scatter(X_TOutx(end),X_TOuty(end),50,'MarkerEdgeColor',[0.1 0.7 0.1],'MarkerFaceColor',[0.1 0.7 0.1],'LineWidth',1)
%     xlabel('East [m]')
%     ylabel('North [m]')
%     grid on
%     legend('Evader Initial Location','Pursuer Initial Location','Evader Trajectory','Pursuer Trajectroy','Catch Location','Location','best')
%     title(sprintf('Problem 1: Evader and Pursuer Trajectories with V_P = %d',V_M))


    % take beta dot
    beta_dotOut = getBetaDot (statesOut,V_T,V_M);

    % take beta double dot
    beta_ddotOut = getBetaDDot (statesOut,V_T,V_M,beta_dotOut);

    if i==1
        statesOut1 = statesOut;
        tOut1 = tOut;
        beta_dotOut1 = beta_dotOut;
        beta_ddotOut1 = beta_ddotOut;
    elseif i==2
        statesOut2 = statesOut;
        tOut2 = tOut;
        beta_dotOut2 = beta_dotOut;
        beta_ddotOut2 = beta_ddotOut;
    elseif i==3
        statesOut3 = statesOut;
        tOut3 = tOut;
        beta_dotOut3 = beta_dotOut;
        beta_ddotOut3 = beta_ddotOut;
    end


end


for j=1:3
    l1 = size(statesOut1,1);
    l2 = size(statesOut2,1);
    l3 = size(statesOut3,1);
    l = [l1;l2;l3];
    [M,I]=max(l);
end
I

figure(i+1)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
plot(statesOut1(:,3), statesOut1(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'Pursuer Trajectroy with V_P = 6 m/s','Pursuer Trajectroy with V_P = 8 m/s','Pursuer Trajectroy with V_P = 11 m/s',...
    'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s',...
    'Location','best');
title('Problem 1: Evader and Pursuer Trajectories with V_P');

figure(i+2)
plot(timeVec(1:length(beta_dotOut1)-1),abs(beta_dotOut1(1:end-1)),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(timeVec(1:length(beta_dotOut2)-1),abs(beta_dotOut2(1:end-1)),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(timeVec(1:length(beta_dotOut3)-1),abs(beta_dotOut3(1:end-1)),'LineWidth',2,'Color',[0.7 0.6 0.3])
xlabel('$time [s]$','interpreter','latex')
ylabel('$\dot \beta [rad/sec]$','interpreter','latex');
grid on
legend('$\dot \beta$ with $V_P$ = 6 [m/s]','$\dot \beta$ with $V_P$ = 8 [m/s]','$\dot \beta$ with $V_P$ = 11 [m/s]','interpreter','latex','Location','best');
title('$\dot \beta$','interpreter','latex');

figure(i+3)
plot(timeVec(1:length(beta_ddotOut1)-1),abs(beta_ddotOut1(1:end-1)),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(timeVec(1:length(beta_ddotOut2)-1),abs(beta_ddotOut2(1:end-1)),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(timeVec(1:length(beta_ddotOut3)-1),abs(beta_ddotOut3(1:end-1)),'LineWidth',2,'Color',[0.7 0.6 0.3])
xlabel('$time [s]$','interpreter','latex')
ylabel('$\ddot \beta [rad/sec^{2}]$','interpreter','latex');
grid on
legend('$\ddot \beta$ with $V_P$ = 6 [m/s]','$\ddot \beta$ with $V_P$ = 8 [m/s]','$\ddot \beta$ with $V_P$ = 11 [m/s]','interpreter','latex','Location','best');
title('$\ddot \beta$','interpreter','latex');


figure(i+4)
plot(tOut1, statesOut1(:,6),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(tOut2, statesOut2(:,6),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(tOut3, statesOut3(:,6),'LineWidth',2,'Color',[0.7 0.6 0.3])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('Range with V = 6','Range with V = 8','Range with V = 11',...
    'Location','best')
title('Problem 1: Ranges')