function [estimate, samples, evaluations] = ISmonteCarlo2(runNumber,numSteps,seed,pThresholdf,pThresholdq,distLimitf,distLimitq,sumType)


fsteps = generateRandomWalkSteps(runNumber,numSteps,pThresholdf,seed);
frandomWalk=generateRandomWalk(fsteps);
%[fevaluation, fevaluationarray]=evaluatorFunc(frandomWalk,distLimitf);
% f_x = prob_x(fsteps,pThresholdf);


qsteps = generateRandomWalkSteps(runNumber,numSteps,pThresholdq,seed);
qrandomWalk=generateRandomWalk(qsteps);
% [qevaluation, qevaluationarray]=evaluatorFunc(proposalFunc,distLimitq);
q_x = prob_x(qsteps,pThresholdq);
% plot(qrandomWalk','LineWidth',1.5,'Color',[0, 0.5, 0.5, 0.2])
% hist(qrandomWalk(:,end))
f_x = prob_x(qsteps,pThresholdf);


% g_x = indFunc(qrandomWalk(:,end),distLimitq);
g_x =evaluatorFunc(qrandomWalk(:,end),distLimitq);

% g_xThreshold = 1;
% [gevaluation, gevaluationarray]=evaluatorFunc(g_x,g_xThreshold);

likelihoodfq   = f_x./q_x;
plikelihoodfq  = prod(likelihoodfq');
evaluations    = g_x.*plikelihoodfq';


samples = qrandomWalk; %???
%evaluations = sum(evaluations)
%evaluations = sum(evaluatorFunc(g_x.*f_x./q_x,distLimitq));


% 
% steps = generateRandomWalkSteps(runNumber,numSteps,pThreshold,seed);
% randomWalk=generateRandomWalk(steps);
% samples = randomWalk;
% evaluations=evaluatorFunc(randomWalk,distLimit);

if sumType == 0
    estimate = sum(evaluations) / (runNumber); % indicator func sum columnwise summation
else
    estimate = cumsum(evaluations) / (1:runNumber+1); % columnwise summation
end
end