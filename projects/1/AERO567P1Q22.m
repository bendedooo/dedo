% 567 3D RW
clc
clear all
close all

seedX = 1333; % rng seed for reproducibility
seedY = 4806; % rng seed for reproducibility
seedZ = 1881; % rng seed for reproducibility
seed  = [seedX seedY seedZ];
% rng seed fix
% s = rng;
% s.Seed = seed;
% rng(s.Seed)

%% PART A

numSteps  = 100;
runNumber = 1000;
stepsX = generateRandomWalkStepsNormal(runNumber,numSteps,seed(1));
stepsY = generateRandomWalkStepsNormal(runNumber,numSteps,seed(2));
stepsZ = generateRandomWalkStepsNormal(runNumber,numSteps,seed(3));
randomWalk = cat(3,stepsX,stepsY,stepsZ);

randomWalkT = cat(3,(randomWalk(:,:,1))',(randomWalk(:,:,2))',(randomWalk(:,:,3))');
posEnd = sum(randomWalkT);
% xyzEnd = squeeze(posEnd(:,end,:))
posEnd = squeeze(posEnd);
distEnd = vecnorm(posEnd'); % If A is a matrix, then vecnorm returns the norm of each column.

% figure(1)
% plot(randomWalk','LineWidth',1.5,'Color',[0, 0.5, 0.5, 0.2])
% title('Question 2.1: Part A','Random Walk')
% xlabel('N')
% ylabel('S')
% legend('Random Walk_n','Location','best')
% grid on

figure(2)
hist(distEnd)
% h.FaceColor = [0 0.5 0.5];
title('Question 2.2: Part A','Random Walk Distance Histogram')
xlabel('S_N')
ylabel('Frequency')
legend('Random Walk_n','Location','best')

%% PART B
distLimit = 10;
runNumber = 100000;
sumType = 0; 

[estimate, samples, evaluations] = monteCarloNormal(runNumber,numSteps,seed,distLimit,sumType);
EvalS10=evaluations;
estimate

% PS10 = estimate;
% EvalS10 = evaluations;

%% PART C
% IS Sampling Rules
distLimitc = 55;
runNumber = 100000;
sumType = 0; 

mu1 = 0;
sigma1 = 1;

mu2 = 0.25;
sigma2 = 1;

[estimate, samples, evaluations] = ISmonteCarloNormal(runNumber,numSteps,seed,mu1,sigma1,mu2,sigma2,distLimitc,sumType);
EvalS55 = evaluations;
estimate

% PS55 = estimate;
% EvalS55 = evaluations;

%% PART D

stdError10 = std(EvalS10)/sqrt(runNumber);

stdError55 = std(EvalS55)/sqrt(runNumber);

z = 1.96;
% prob 0.95 ==> z score 1.96

meanEvalS10 = mean(EvalS10);
varEvalS10 = var(EvalS10);
intervalS10 = [meanEvalS10-z*varEvalS10/(sqrt(runNumber)) meanEvalS10-z*varEvalS10/(sqrt(runNumber))];

       
meanEvalS55 = mean(EvalS55);
varEvalS55 = var(EvalS55);
intervalS55 = [meanEvalS55-z*varEvalS55/(sqrt(runNumber)) meanEvalS55-z*varEvalS55/(sqrt(runNumber))];



%% extra
AnalyticalP10 = 1 - gamcdf(distLimit,3/2,6)
AnalyticalP55 = 1 - gamcdf(distLimitc,3/2,6)


%[estimate, samples, evaluations] = monteCarloNormalIS(runNumber,numSteps,seed,distLimit,sumType);
%evaluations

