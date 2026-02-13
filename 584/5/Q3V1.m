% HW 5 P2: Constant Bearig Puruit and Proportional
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
V_M = 6;
% theta_T = pi/2; % bunu formulize et

lambdaArray= [0.25 0.5 0.75 0.9 1, 2, 5, 50]; % for 8 (gamma>3/2), for 11 (gamma>2)

R0 = sqrt((X_T0x-X_M0x)^2+(X_T0y-X_M0y)^2); % Range
beta0 = atan2((X_T0y-X_M0y),(X_T0x-X_M0x));
theta0 = beta0;
timeVec = 0:0.01:100;
states0 = [X_M0; X_T0; beta0; R0];


odeFunc = @(t,x)CBP(t,x,V_M,V_T);
optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
[tOutCBP, statesOutCBP] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);

states0 = [X_M0; X_T0; beta0; R0; theta0];


statesOutMatrix = [];

for i = 1:length(lambdaArray)
    lambda=lambdaArray(i);
    odeFunc = @(t,x)PP(t,x,V_M,V_T,lambda);

    optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
    [tOut, statesOut] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);


    if i==1
        statesOut1 = statesOut;
        tOut1 = tOut;
    elseif i==2
        statesOut2 = statesOut;
        tOut2 = tOut;
    elseif i==3
        statesOut3 = statesOut;
        tOut3 = tOut;
    elseif i==4
        statesOut4 = statesOut;
        tOut4 = tOut;
    elseif i==5
        statesOut5 = statesOut;
        tOut5 = tOut;
    elseif i==6
        statesOut6 = statesOut;
        tOut6 = tOut;
    elseif i==7
        statesOut7 = statesOut;
        tOut7 = tOut;
    elseif i==8
        statesOut8 = statesOut;
        tOut8 = tOut;

    end


end


% figure(i+1)
% scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
% hold on
% scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% plot(statesOut1(:,3), statesOut1(:,4),'LineWidth',2,'Color',[.5 0 .5])
% plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
% plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
% plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
% plot(statesOut4(:,1), statesOut4(:,2),'LineWidth',2,'Color',[0.7 0.6 0.1])
% % scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% % scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
% % scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
% % scatter(statesOut4(end,1),statesOut4(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.1],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
% xlabel('East [m]')
% ylabel('North [m]')
% grid on
% legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
%     'Pursuer Trajectroy with V_P = 6 m/s','Pursuer Trajectroy with V_P = 8 m/s','Pursuer Trajectroy with V_P = 11 m/s','Pursuer Trajectroy with V_P = 11 m/s',...
%     'Location','best')
% % 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
% title(sprintf('Problem 1: Evader and Pursuer Trajectories with V_P = %d',V_M))

%% plots
figure(1)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
plot(statesOut1(:,3), statesOut1(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(statesOut4(:,1), statesOut4(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
scatter(statesOut4(end,1),statesOut4(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'Pursuer Trajectroy with \lambda = 0.25','Pursuer Trajectroy with \lambda = 0.50','Pursuer Trajectroy with \lambda = 0.75','Pursuer Trajectroy with \lambda = 0.9',...
    'Location','best')
% 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
title('Problem 3: PP Evader and Pursuer Trajectories with \lambda={0.25,0.5,0.75,0.9}')

figure(2)
plot(tOut1(:,1), statesOut1(:,6),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(tOut2(:,1), statesOut2(:,6),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(tOut3(:,1), statesOut3(:,6),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(tOut4(:,1), statesOut4(:,6),'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('Range with \lambda = 0.25','Range with \lambda = 0.50','Range with \lambda = 0.75','Range with \lambda = 0.9',...
    'Location','best')
% 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
title('Problem 3: PP Range with \lambda={0.25,0.5,0.75,0.9}')

figure(3)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
plot(statesOut5(:,3), statesOut5(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut5(:,1), statesOut5(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
plot(statesOut6(:,1), statesOut6(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(statesOut7(:,1), statesOut7(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(statesOut8(:,1), statesOut8(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
plot(statesOutCBP(:,1), statesOutCBP(:,2),'LineWidth',2,'Color',[0.1 0.1 0.1])
scatter(statesOut5(end,1),statesOut5(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
scatter(statesOut6(end,1),statesOut6(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
scatter(statesOut7(end,1),statesOut7(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
scatter(statesOut8(end,1),statesOut8(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
scatter(statesOutCBP(end,1),statesOutCBP(end,2),50,'MarkerEdgeColor',[0.1 0.1 0.1],'MarkerFaceColor',[0.1 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'PP Pursuer Trajectroy with \lambda = 1','PP Pursuer Trajectroy with \lambda = 2','PP Pursuer Trajectroy with \lambda = 5','PP Pursuer Trajectroy with \lambda = 50','CBP Pursuer Trajectroy',...
    'Location','best')
% 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
title('Problem 3: CBP and PP Evader and Pursuer Trajectories with \lambda={1,2,5,50}')

figure(4)
plot(tOut5, statesOut5(:,6),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(tOut6, statesOut6(:,6),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(tOut7, statesOut7(:,6),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(tOut8, statesOut8(:,6),'LineWidth',2,'Color',[0.7 0.1 0.1])
plot(tOutCBP, statesOutCBP(:,6),'LineWidth',2,'Color',[0.1 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('Range with \lambda = 1','Range with \lambda = 2','Range with \lambda = 5','Range with \lambda = 50','CBP Range',...
    'Location','best')
title('Problem 3: CBP and PP Range with \lambda={1,2,5,50}')