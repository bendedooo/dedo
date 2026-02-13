% 4.1 deneme

clear all
close all
clc

seed = 1333; % rng seed for reproducibility
rng; %seed fix
s = rng;
s.Seed = seed;
rng(s.Seed);

mu = 0.05;
sigma = 0.2;
Y_0 = 1;
% 0<t<1
tFinal = 1;
% dt =1/N;
dt = 0.05;
% N = 1/dt;
runNumber = 100;


%% PART 1

W = bms(tFinal, dt, runNumber);
Y = gbms(W,mu,sigma,Y_0,dt);

% plots

timeVec = 0:dt:tFinal;

figure(1)
plot(timeVec,W')
title('Question 4.1: Part 1','Wiener Process ')
xlabel('Time')
ylabel('Walked Distance')
legend('Wiener Process_n','Location','best')
grid on

figure(2)
plot(timeVec,Y')
title('Question 4.1: Part 1','Geometric Brownian Motion')
xlabel('Time')
ylabel('Walked Distance')
legend('Geometric Brownian Motion_n','Location','best')
grid on


var_Y=var(Y(:,end));
exp_Y=mean(Y(:,end));


AnalyticExp_Y = exp(mu - (sigma^2)/2 + (sigma^2)/2 );
AnalyticVar_Y = exp(2*mu - sigma^2 + 2*sigma^2)- exp(2*mu - sigma^2 + sigma^2);

%% PART 2
% clear all
% 
% seed = 1333; % rng seed for reproducibility
% rng %seed fix
% s = rng;
% s.Seed = seed;
% rng(s.Seed);

% tFinal = 1;
% dt = 0.05;     % fine time step

runNumber = 1;
m=4;
dtm = dt*m;    % coarse time step


Wf = bms(tFinal, dt, runNumber);
Wc = bmf2c(tFinal, dtm, runNumber, Wf, m);

Yf = gbms(Wf,mu,sigma,Y_0,dt);
Yc = gbms(Wc,mu,sigma,Y_0,dtm);

% plots

timeVecC = 0:dtm:tFinal;
timeVecF = 0:dt:tFinal;

figure(3)
plot(timeVecF,Wf,'m-x','LineWidth',2,'MarkerSize',8)
hold on
plot(timeVecC,Wc,'b-o','LineWidth',2,'MarkerSize',8)
title('Question 4.1: Part 2','Wiener Process vs dt')
xlabel('Time')
ylabel('Walked Distance')
legend('Fine Wiener Process','Coarse Wiener Process','Location','best')
grid on

figure(4)
plot(timeVecF,Yf,'m-x','LineWidth',2,'MarkerSize',8)
hold on
plot(timeVecC,Yc,'b-o','LineWidth',2,'MarkerSize',8)
title('Question 4.1: Part 2','Geometric Brownian Motion vs dt')
xlabel('Time')
ylabel('Walked Distance')
legend('Fine Geometric Brownian Motion','Coarse Geometric Brownian Motion','Location','best')
grid on

%% PART 3

runNumber = 1000;
i = [2 3 4 5];

dt = 4^-i(4);
Wf = bms(tFinal, dt, runNumber);
Yf = gbms(Wf,mu,sigma,Y_0,dt);

m1=4;
dtc1 = dt*m1;    % coarse time step 1
Wc1 = bmf2c(tFinal, dtc1, runNumber, Wf, m1);
Yc1 = gbms(Wc1,mu,sigma,Y_0,dtc1);

m2=4^2;
dtc2 = dt*m2;
Wc2 = bmf2c(tFinal, dtc2, runNumber, Wf, m2);
Yc2 = gbms(Wc2,mu,sigma,Y_0,dtc2);

m3=4^3;
dtc3 = dt*m3;
Wc3 = bmf2c(tFinal, dtc3, runNumber, Wf, m3);
Yc3 = gbms(Wc3,mu,sigma,Y_0,dtc3);

% plots

timeVecF = 0:dt:tFinal;
timeVecC1 = 0:dtc1:tFinal;
timeVecC2 = 0:dtc2:tFinal;
timeVecC3 = 0:dtc3:tFinal;

% figure(5)
% plot(timeVecF,Wf,'m-x','LineWidth',2,'MarkerSize',8)
% hold on
% plot(timeVecC1,Wc1,'b-o','LineWidth',2,'MarkerSize',8)
% plot(timeVecC2,Wc2,'g-o','LineWidth',2,'MarkerSize',8)
% plot(timeVecC2,Wc2,'r-o','LineWidth',2,'MarkerSize',8)
% title('Question 4.1: Part 2','Wiener Process vs dt')
% xlabel('Time')
% ylabel('Walked Distance')
% legend('Fine Wiener Process','Coarse Wiener Process #1','Coarse Wiener Process #2','Coarse Wiener Process #3','Location','best')
% grid on
% 
% figure(6)
% plot(timeVecF,Yf,'m-x','LineWidth',2,'MarkerSize',8)
% hold on
% plot(timeVecC1,Yc1,'b-o','LineWidth',2,'MarkerSize',8)
% plot(timeVecC2,Yc2,'g-o','LineWidth',2,'MarkerSize',8)
% plot(timeVecC2,Yc2,'r-o','LineWidth',2,'MarkerSize',8)
% title('Question 4.1: Part 2','Geometric Brownian Motion vs dt')
% xlabel('Time')
% ylabel('Walked Distance')
% legend('Fine Geometric Brownian Motion','Coarse Geometric Brownian Motion #1','Coarse Geometric Brownian Motion #2','Coarse Geometric Brownian Motion #3','Location','best')
% grid on

Yf_end = Yf(:,end);
Yc1_end = Yc1(:,end);
Yc2_end = Yc2(:,end);
Yc3_end = Yc3(:,end);

% E[Yf] = E[Yc1]+E[Yf-Yc1]
% E[Yf] = E[Yc2]+E[Yc1-Yc2]+E[Yf-Yc1]
% E[Yf] = E[Yc3]+E[Yc2-Yc3]+E[Yc1-Yc2]+E[Yf-Yc1]



V_Yc3 = var(Yc3_end);
V_Yc2Yc3 = var(Yc2_end-Yc3_end);
V_Yc1Yc2 = var(Yc1_end-Yc2_end);
V_YfYc1 = var(Yf_end-Yc1_end);

VarY = [V_Yc3 V_Yc2Yc3 V_Yc1Yc2 V_YfYc1];

V_Yf = V_Yc3 + V_Yc2Yc3 + V_Yc1Yc2 + V_YfYc1;

E_Yc3    = sum(Yc3_end)/length(Yc3_end);
E_Yc2Yc3 = sum(Yc2_end-Yc3_end)/length(Yc2_end-Yc3_end);
E_Yc1Yc2 = sum(Yc1_end-Yc2_end)/length(Yc1_end-Yc2_end);
E_YfYc1  = sum(Yf_end-Yc1_end)/length(Yf_end-Yc1_end);

E_Yf = E_Yc3 + E_Yc2Yc3 + E_Yc1Yc2 + E_YfYc1;
EY   = [E_Yc3  E_Yc2Yc3  E_Yc1Yc2  E_YfYc1];

% P 

P_Yc3 = exp(-0.05)*max(0,(Yc3_end)-1);
P_Yc2 = exp(-0.05)*max(0,(Yc2_end)-1);
P_Yc1 = exp(-0.05)*max(0,(Yc1_end)-1);
P_Yf = exp(-0.05)*max(0,(Yf_end)-1);

P_Yc2Yc3 = P_Yc2-P_Yc3;
P_Yc1Yc2 = P_Yc1-P_Yc2;
P_YfYc1  = P_Yf-P_Yc1;
% 
% P_Yc2Yc3 = exp(-0.05)*max(0,(Yc2_end-Yc3_end)-1);
% P_Yc1Yc2 = exp(-0.05)*max(0,(Yc1_end-Yc2_end)-1);
% P_YfYc1 = exp(-0.05)*max(0,(Yf_end-Yc1_end)-1);

% P Expectation

EP_Yc3    = sum(P_Yc3)/length(P_Yc3);
EP_Yc2Yc3 = sum(P_Yc2Yc3)/length(P_Yc2Yc3);
EP_Yc1Yc2 = sum(P_Yc1Yc2)/length(P_Yc1Yc2);
EP_YfYc1  = sum(P_YfYc1)/length(P_YfYc1);

EP_Yf = EP_Yc3 + EP_Yc2Yc3 + EP_Yc1Yc2 + EP_YfYc1;
EYP   = [EP_Yc3 EP_Yc2Yc3 EP_Yc1Yc2 EP_YfYc1];
% P Variance

VP_Yc3 = var(P_Yc3);
VP_Yc2Yc3 = var(P_Yc2Yc3);
VP_Yc1Yc2 = var(P_Yc1Yc2);
VP_YfYc1  = var(P_YfYc1);

VP_Yf = VP_Yc3 + VP_Yc2Yc3 + VP_Yc1Yc2 + VP_YfYc1;
VarYP = [VP_Yc3 VP_Yc2Yc3 VP_Yc1Yc2 VP_YfYc1];

figure(7)
bar((VarY))
set(gca,'YScale','log')
title('Question 4.1: Part 3','Geometric Brownian Motion MLMC Variance vs Levels')
xlabel('Levels')
ylabel('log(Variance)')
legend('Geometric Brownian Motion Variances','Coarse Geometric Brownian Motion #1','Coarse Geometric Brownian Motion #2','Coarse Geometric Brownian Motion #3','Location','best')
grid on


figure(8)
bar((VarYP))
set(gca,'YScale','log')
title('Question 4.1: Part 3','Geometric Brownian Motion Pay Off Function MLMC Variance vs Levels')
xlabel('Levels')
ylabel('log(Variance)')
legend('Geometric Brownian Motion Variances','Coarse Geometric Brownian Motion #1','Coarse Geometric Brownian Motion #2','Coarse Geometric Brownian Motion #3','Location','best')
grid on



figure(9)
bar((EY))
set(gca,'YScale','log')
title('Question 4.1: Part 3','Geometric Brownian Motion MLMC Expectation vs Levels')
xlabel('Levels')
ylabel('log(Expectation)')
legend('Geometric Brownian Motion Estimates','Coarse Geometric Brownian Motion #1','Coarse Geometric Brownian Motion #2','Coarse Geometric Brownian Motion #3','Location','best')
grid on


figure(10)
bar(((EYP)))
set(gca,'YScale','log')
title('Question 4.1: Part 3','Geometric Brownian Motion Pay Off Function MLMC Expectation vs Levels')
xlabel('Levels')
ylabel('log(Expectation)')
% ylim([10^-6 10^-1])
legend('Geometric Brownian Motion Estimates','Coarse Geometric Brownian Motion #1','Coarse Geometric Brownian Motion #2','Coarse Geometric Brownian Motion #3','Location','best')
grid on

%% PART 4

t4=0.029989;
t3= 0.032751;
t2=0.041187;
t1=0.054871;

t= t1+t2+t3+t4;

normalized_t4 = t4/t;
normalized_t3 = t3/t;
normalized_t2 = t2/t;
normalized_t1 = t1/t;

tarVar=0.005;
lambda = tarVar^(-2)*(sqrt(normalized_t4*VP_Yc3)*sqrt(normalized_t3*VP_Yc2Yc3)*sqrt(normalized_t2*VP_Yc1Yc2)*sqrt(normalized_t1*VP_YfYc1))^2;

N_L = lambda.*sqrt([VP_Yc3 VP_Yc2Yc3 VP_Yc1Yc2 VP_YfYc1].\[normalized_t4 normalized_t3 normalized_t2 normalized_t1]);


%% cop
% % fine
% 
% mf = 1;
% Wf = bms(tFinal, mf*dt, runNumber);
% 
% Wf=Wf';
% sizeWf = size(Wf);
% Yf(1:sizeWf(1),1) = Y_0;
% 
% for i=1:sizeWf(2)-1
%     b = mu * Yf(:,i);
%     h = sigma * Yf(:,i);
%     dWf = Wf(:,i+1)-Wf(:,i);
%     Yf(:,i+1) = Yf(:,i) + b.*mf*dt + h.*dWf;
% end
% plot(Wf')
% figure(1)
% plot(Yf')
% 
% % coarse
% Wc = Wf(1:4:end,:)
% % mc = 4;
% % Wc = bms(tFinal, mc*dt, runNumber);
% % Wc=Wc';
% sizeWc = size(Wc);
% Yc(1:sizeWc(1),1) = Y_0;
% 
% for i=1:sizeWc(2)-1
%     b = mu * Yc(:,i);
%     h = sigma * Yc(:,i);
%     dWc = Wc(:,i+1)-Wc(:,i);
%     Yc(:,i+1) = Yc(:,i) + b.*mc*dt + h.*dWc;
% end
% 
% plot(Wc')
% figure(2)
% plot(Yc')
% 
% 
% figure(3)
% plot(Yc(1,:))
% hold on 
% plot(Yf(1,:))
% 
% 
% 
% % dw = w_tn+1 - w_tn
% % Y_n+1 = Y_n + b * dt + h * dw
% 
% % pay off function P=exp(-0.05)max(0,Y(1)-1)

%% dersteki ornek

% a=bms(1,0.01,500)
% %samples = [0  0]
% plot(a)



% dw = w_tn+1 - w_tn
% Y_n+1 = Y_n + b * dt + h * dw

% pay off function P=exp(-0.05)max(0,Y(1)-1)
