% HW 5 P2: Constant Bearig Pursuit and Direct Pursuit
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

    odeFunc = @(t,x)CBP(t,x,V_M,V_T);

    optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
    [tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);

    odeFuncDP = @(t,x)Rdotbetadot(t,x,V_M,V_T);


    optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
    [tOutDP, statesOutDP] = ode45(@(t,x)odeFuncDP(t,x),timeVec, states0,optPos);


    %     X_TOutx = statesOut(:,3);
    %     X_TOuty = statesOut(:,4);
    %     X_MOutx = statesOut(:,1);
    %     X_MOuty = statesOut(:,2);
    %
    %
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

    if i==1
        statesOut1 = statesOut;
        statesOutDP1 = statesOutDP;
        tOut1 = tOut;
        tOutDP1 = tOutDP;
    elseif i==2
        statesOut2 = statesOut;
        statesOutDP2 = statesOutDP;
        tOut2 = tOut;
        tOutDP2 = tOutDP;

    elseif i==3
        statesOut3 = statesOut;
        statesOutDP3 = statesOutDP;
        tOut3 = tOut;
        tOutDP3 = tOutDP;

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

% figure(i+1)
% scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
% hold on
% scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% plot(statesOut1(:,3), statesOut1(:,4),'LineWidth',2,'Color',[.5 0 .5])
% plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
% plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
% plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
% scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
% scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
% xlabel('East [m]')
% ylabel('North [m]')
% grid on
% legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
%     'Pursuer Trajectroy with V_P = 6 m/s','Pursuer Trajectroy with V_P = 8 m/s','Pursuer Trajectroy with V_P = 11 m/s',...
%     'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s',...
%     'Location','best')
% title('Problem 2: Evader and Pursuer Trajectories with V_P')
%% plots
figure(1)
subplot(3,1,1)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
plot(statesOutDP1(:,3), statesOutDP1(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
plot(statesOutDP1(:,1), statesOutDP1(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
scatter(statesOutDP1(end,1),statesOutDP1(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'CBP Pursuer Trajectroy with V_P = 6 m/s','DP Pursuer Trajectroy with V_P = 6 m/s',...
    'Location','best')
title('Problem 2: Evader and Pursuer Trajectories with V_P')
%     'Catch Location with V_P = 6 m/s','DP Catch Location with V_P = 6 m/s',...


subplot(3,1,2)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
plot(statesOutDP2(:,3), statesOutDP2(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(statesOutDP2(:,1), statesOutDP2(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
scatter(statesOutDP2(end,1),statesOutDP2(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'CBP Pursuer Trajectroy with V_P = 8 m/s','DP Pursuer Trajectroy with V_P = 8 m/s',...
    'Location','best')


subplot(3,1,3)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
plot(statesOutDP3(:,3), statesOutDP3(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(statesOutDP3(:,1), statesOutDP3(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
scatter(statesOutDP3(end,1),statesOutDP3(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'CBP Pursuer Trajectroy with V_P = 11 m/s','DP Pursuer Trajectroy with V_P = 11 m/s',...
    'Location','best')

figure(2)
subplot(3,1,1)
plot(tOut1,statesOut1(:,6),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(tOutDP1,statesOutDP1(:,6),'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('CBP Range with V_P = 6 m/s','DP Range with V_P = 6 m/s',...
    'Location','best')
title('Problem 2: Range with V_P')
subplot(3,1,2)
plot(tOut2,statesOut2(:,6),'LineWidth',2,'Color',[0.1 0.6 0.3])
hold on
plot(tOutDP2,statesOutDP2(:,6),'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('CBP Range with V_P = 8 m/s','DP Range with V_P = 8 m/s',...
    'Location','best')
subplot(3,1,3)
plot(tOut3,statesOut3(:,6), 'LineWidth',2,'Color',[0.7 0.6 0.3])
hold on
plot(tOutDP3,statesOutDP3(:,6), 'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('CBP Range with V_P = 11 m/s','DP Range with V_P = 11 m/s',...
    'Location','best')



%% gereksiz plot

% plotlari ayarla