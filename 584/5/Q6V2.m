% HW 5 P6: Noisy DP and PP
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

lambdaArray= [1,1.5, 2, 2.5]; % for 8 (gamma>3/2), for 11 (gamma>2)

R0 = sqrt((X_T0x-X_M0x)^2+(X_T0y-X_M0y)^2); % Range
beta0 = atan2((X_T0y-X_M0y),(X_T0x-X_M0x));
timeVec = 0:0.01:100;

states0 = [X_M0; X_T0; beta0; R0];
% odeFunc = @(t,x)DP(t,x,V_M,V_T);
optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
% [tOutDP, statesOutDP] = ode45(@(t,x)odeFunc(t,x),timeVec, states0,optPos);


% theta0 = beta0;
% states0 = [X_M0; X_T0; beta0; R0; theta0];



for i = 1:length(lambdaArray)
    lambda=lambdaArray(i);
    %     states = states0';
    if i ==1
        theta0=0;
    else
        theta0 = beta0;
    end
    states0 = [X_M0; X_T0; beta0; R0; theta0];
    states = states0';
    statesOut = zeros(length(timeVec),length(states0));
    statesOut(1,:)=states0';


    for ii = 2:length(timeVec)
        
        beta = states(end,5);
        betaMeas = beta+(0.25)*randn;
        theta  = lambda*betaMeas+theta0;
        odeFunc = @(t,x)PP(t,x,V_M,V_T,lambda,theta);

        optPos = odeset('Events', @StopEventPos,'RelTol', 1e-6, 'AbsTol', 1e-6);
        [tout, states,~,~,ie] = ode45(@(t,x)odeFunc(t,x),[ii-2 ii-1]*0.01, states(end,:),optPos);
        statesOut(ii,:) = states(end,:)';
        if ~isempty(ie)
            break
        end
        
    end

        tOut = 0:0.01:0.01*(ii-1);%timeVec(1:ii+1);
        statesOut(ii+1:end,:)=[];


    if i==1
        statesOutDP = statesOut;
        tOutDP = tOut;
    elseif i==2
        statesOut1 = statesOut;
        tOut1 = tOut;
    elseif i==3
        statesOut2 = statesOut;
        tOut2 = tOut;
    elseif i==4
        statesOut3 = statesOut;
        tOut3 = tOut;

    end


end

%% plots

figure(1)
scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
plot(statesOut3(:,3), statesOut3(:,4),'LineWidth',2,'Color',[.5 0 .5])
plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(statesOutDP(:,1), statesOutDP(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
scatter(statesOutDP(end,1),statesOutDP(end,2),50,'MarkerEdgeColor',[0.7 0.1 0.1],'MarkerFaceColor',[0.7 0.1 0.1],'LineWidth',1)
xlabel('East [m]')
ylabel('North [m]')
grid on
legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
    'PP Pursuer Trajectroy with \lambda = 1.5','PP Pursuer Trajectroy with \lambda = 2','PP Pursuer Trajectroy with \lambda = 2.5','DP Pursuer Trajectroy',...
    'Location','best')
% 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
title('Problem 6: DP and PP Evader and Pursuer Trajectories with Noisy Measurements')

figure(2)
plot(tOut1', statesOut1(:,6),'LineWidth',2,'Color',[0, 0.5, 0.5])
hold on
plot(tOut2', statesOut2(:,6),'LineWidth',2,'Color',[0.1 0.6 0.3])
plot(tOut3', statesOut3(:,6),'LineWidth',2,'Color',[0.7 0.6 0.3])
plot(tOutDP', statesOutDP(:,6),'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time [s]')
ylabel('Range [m]')
grid on
legend('Range with \lambda = 1.5','Range with \lambda = 2','Range with \lambda = 2.5','CBP Range',...
    'Location','best')
% 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
title('Problem 6: DP and PP Range with \lambda={1.5,2,2.5}')

%% plots cop
% 
% figure(i+1)
% scatter(X_T0x,X_T0y,50,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
% hold on
% scatter(X_M0x,X_M0y,50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% plot(statesOut1(:,3), statesOut1(:,4),'LineWidth',2,'Color',[.5 0 .5])
% plot(statesOut1(:,1), statesOut1(:,2),'LineWidth',2,'Color',[0, 0.5, 0.5])
% plot(statesOut2(:,1), statesOut2(:,2),'LineWidth',2,'Color',[0.1 0.6 0.3])
% plot(statesOut3(:,1), statesOut3(:,2),'LineWidth',2,'Color',[0.7 0.6 0.3])
% plot(statesOutDP(:,1), statesOutDP(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
% % plot(statesOut4(:,1), statesOut4(:,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
% % scatter(statesOut1(end,1),statesOut1(end,2),50,'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0, 0.5, 0.5],'LineWidth',1)
% % scatter(statesOut2(end,1),statesOut2(end,2),50,'MarkerEdgeColor',[0.1 0.6 0.3],'MarkerFaceColor',[0.1 0.6 0.3],'LineWidth',1)
% % scatter(statesOut3(end,1),statesOut3(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.3],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
% % scatter(statesOut4(end,1),statesOut4(end,2),50,'MarkerEdgeColor',[0.7 0.6 0.1],'MarkerFaceColor',[0.7 0.6 0.3],'LineWidth',1)
% xlabel('East [m]')
% ylabel('North [m]')
% grid on
% legend('Evader Initial Location','Pursuer Initial Location','Northernmost Evader Trajectory',...
%     'PP Pursuer Trajectroy with \lambda = 1.5','PP Pursuer Trajectroy with \lambda = 2','PP Pursuer Trajectroy with \lambda = 2.5','DP Pursuer Trajectroy',...
%     'Location','best')
% % 'Catch Location with V_P = 6 m/s','Catch Location with V_P = 8 m/s','Catch Location with V_P = 11 m/s','Catch Location with V_P = 11 m/s',...
% title(sprintf('Problem 1: Evader and Pursuer Trajectories with V_P = %d',V_M))
