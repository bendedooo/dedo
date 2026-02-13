% HW 4 P3
clc
clear all
close all
format long

load('rcwA.mat');
actualTrajectory = rcwA_Ts_0_01;

g = 9.80665;
gVec = [0 0 -g]'; % in frame A

k= 1:1:2000;
T = 0.01;
phi_0 = pi/6;


D1 = 0.1*eye(3);
D2 = 0.1*eye(3);


omega = [0; 0; 1];



att_0 = [phi_0 0 0]; % [phi; theta; psi]
C_A2B_0 = Angle2DCM_1(phi_0);
v_CW_0 = rdot(phi_0)';
r_CW_0 = [1 0 0]';
phi = phi_0;
% states [r rdot], [x y z]
state0 = [r_CW_0; v_CW_0];

%% PART A

% variables to be saved
states =[];
att =[];


% KF Param Init
P00    = rand*eye(6);
P11    = P00;
X11    = state0;
C_A2B1 = C_A2B_0;
% variables to be saved
statesi =[X11];
%deneme
angi = [];
acc =[];
for i=1:length(k) % for all steps
    ki = k(i); % step
    t  = ki*T; % time
    % create acc and gyro measurements
    accMeasOut = accMeas(g,phi,t);
    noise2     = diag(randn(3)).*eye(3);
    accNoisy   = accMeasOut + diag(D2*noise2);
    noise1     = diag(randn(3)).*eye(3);
    omegaNoisy = omega + diag(D2*noise1);

    measAvailable = 0;

    %     m = [m measAvailable];
    [X22, P22, C_A2B2] = KFP3B(X11,C_A2B1,P11,accNoisy, omegaNoisy, zeros(3,1),measAvailable,T);
    % Update KF inits
    X11    = X22;
    P11    = P22;
    C_A2B1 = C_A2B2;
    %     ang=DCM2Angle_321(C_A2B1);
    [yaw, pitch, roll]=dcm2angle(C_A2B1);
    ang=[yaw, pitch, roll]';
    % save variables
    statesi =[statesi X22];
    %deneme
    angi = [angi ang];
    acc = [acc accNoisy];
end

posnA = statesi(1:3,:);

figure(1)
% scatter3(X_MOCAP1(1,2:end), X_MOCAP1(2,2:end), X_MOCAP1(3,2:end),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1)
% hold on
scatter3(posnA(1,1), posnA(2,1), posnA(3,1),100,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
plot3(posnA(1,:), posnA(2,:), posnA(3,:),'LineWidth',1.5,'Color',[0, 0.5, 0.5])
plot3(actualTrajectory(1,:), actualTrajectory(2,:), actualTrajectory(3,:),'LineWidth',2,'Color',[0.7 0.7 0.1])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Quadcopter Trajectory','Quadcopter Initial Location','Actual Trajectory','Location','best')
title('Problem 3: 3D Trajectory of the Quadcopter with Inertial Navigation KF')


%% Part B
% MOCAP Measurements
posMeas = actualTrajectory;
D3 = 0.005*eye(3);
noise3 = diag(randn(3)).*eye(3);
posMeasNoisy = posMeas + diag(D3*noise3);
posMeasNoisy(:,1)=[];

% MOCAP period
T_MOCAP = [1 0.1];
m=[];

% variables to be saved
states =[];
att =[];

for j=1:length(T_MOCAP) % for all MOCAP periods
    % KF Param Init
    P00    = 10*eye(6);
    P11    = P00;
    X11    = state0;
    C_A2B1 = C_A2B_0;
    % variables to be saved
    statesi =[X11];
    %deneme
    angi = [];
    acc =[];
    for i=1:length(k) % for all steps
        ki = k(i); % step
        t  = ki*T; % time
        % create acc and gyro measurements
        accMeasOut = accMeas(g,phi,t);
        noise2     = diag(randn(3)).*eye(3);
        accNoisy   = accMeasOut + diag(D2*noise2);
        noise1     = diag(randn(3)).*eye(3);
        omegaNoisy = omega + diag(D2*noise1);
        % check if MOCAP is available
        if mod(t,T_MOCAP(j))==0
            measAvailable = 1;
        else
            measAvailable = 0;
        end
        %     m = [m measAvailable];
        [X22, P22, C_A2B2] = KFP3B(X11,C_A2B1,P11,accNoisy, omegaNoisy, posMeasNoisy(:,i),measAvailable,T);
        % Update KF inits
        X11    = X22;
        P11    = P22;
        C_A2B1 = C_A2B2;
        %     ang=DCM2Angle_321(C_A2B1);
        [yaw, pitch, roll]=dcm2angle(C_A2B1);
        ang=[yaw, pitch, roll]';
        % save variables
        statesi =[statesi X22];
        %deneme
        angi = [angi ang];
        acc = [acc accNoisy];
    end
    % save variables
    states =cat(3,states,statesi);
    att = cat(3,att,angi);
end

X_MOCAP1=states(1:3,:,1);
X_MOCAP2=states(1:3,:,2);

att_MOCAP1=rad2deg(att(:,:,1)); % bu nasil imaginary gelebilir ya????
att_MOCAP2=rad2deg(att(:,:,2));

figure(2)
% scatter3(X_MOCAP1(1,2:end), X_MOCAP1(2,2:end), X_MOCAP1(3,2:end),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1)
scatter3(X_MOCAP1(1,1), X_MOCAP1(2,1), X_MOCAP1(3,1),100,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
plot3(X_MOCAP1(1,:), X_MOCAP1(2,:), X_MOCAP1(3,:),'LineWidth',1.5,'Color',[0, 0.5, 0.5])
plot3(actualTrajectory(1,:), actualTrajectory(2,:), actualTrajectory(3,:),'LineWidth',2,'Color','r')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Quadcopter Initial Location','Quadcopter Trajectory','Actual Trajectory','Location','best')
title('Problem 3: 3D Trajectory of the Quadcopter with T_{MOCAP}=1 s')

figure(3)
% scatter3(X_MOCAP2(1,2:end), X_MOCAP2(2,2:end), X_MOCAP2(3,2:end),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1)
scatter3(X_MOCAP2(1,1), X_MOCAP2(2,1), X_MOCAP2(3,1),100,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
plot3(X_MOCAP2(1,:), X_MOCAP2(2,:), X_MOCAP2(3,:),'LineWidth',1.5,'Color',[0.1 0.7 0.1])
plot3(actualTrajectory(1,:), actualTrajectory(2,:), actualTrajectory(3,:),'LineWidth',1,'Color','r')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Quadcopter Initial Location','Quadcopter Trajectory','Actual Trajectory','Location','best')
title('Problem 3: 3D Trajectory of the Quadcopter with T_{MOCAP}=0.1 s')


%% plots 

figure(4)
scatter3(X_MOCAP2(1,1), X_MOCAP2(2,1), X_MOCAP2(3,1),100,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
hold on
plot3(X_MOCAP1(1,:), X_MOCAP1(2,:), X_MOCAP1(3,:),'LineWidth',1.5,'Color',[0, 0.5, 0.5])
plot3(X_MOCAP2(1,:), X_MOCAP2(2,:), X_MOCAP2(3,:),'LineWidth',1.5,'Color',[0.1 0.7 0.1])
hold on
plot3(actualTrajectory(1,:), actualTrajectory(2,:), actualTrajectory(3,:),'LineWidth',1,'Color','r');%[0.7 0.7 0.1]);
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
grid on
legend('Quadcopter Initial Location','Quadcopter Trajectory T_{MOCAP}=1 s','Quadcopter Trajectory T_{MOCAP}=0.1 s','Actual Trajectory','Location','best')
title('Problem 3: Part B: 3D Trajectories of the Quadcopter for T_{MOCAP}')



coord =['x','y','z'];
figure(5)
for ir = 1:3
    subplot(3,1,ir);
    plot([0 k*T], X_MOCAP1(ir,:),'LineWidth',2,'Color',[0, 0.5, 0.5]);
    hold on
    plot([0 k*T], X_MOCAP2(ir,:),'LineWidth',2,'Color',[0.1 0.7 0.1]);
    plot([0 k*T], actualTrajectory(ir,:),'LineStyle','-','LineWidth',0.75,'Color','r');%[0.7 0.7 0.1]);
    grid on
    xlabel('time [s]')
    ylabel(sprintf('{%s} coordinate [m]',coord(ir)))
end
legend('Quadcopter Trajectory T_{MOCAP}=1 s','Quadcopter Trajectory T_{MOCAP}=0.1 s','Actual Trajectory','Location','best');
sgtitle('Problem 3: Part B: Components of r_{{c/w}|A} vs Time')




















% 
% figure(3)
% scatter3(actualTrajectory(1,2:end), actualTrajectory(2,2:end), actualTrajectory(3,2:end),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1)
% hold on
% scatter3(actualTrajectory(1,1), actualTrajectory(2,1), actualTrajectory(3,1),100,'MarkerEdgeColor',[.5 0 .5],'MarkerFaceColor',[.5 0 .5],'LineWidth',1)
% %plot3(rOut(1,:), rOut(2,:), rOut(3,:),'LineWidth',1.5,'Color',[0, 0.5, 0.5])
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% legend('Quadcopter Trajectory','Quadcopter Initial Location','Location','best')
% title('Problem 4: 3D Trajectory of the Quadcopter')
