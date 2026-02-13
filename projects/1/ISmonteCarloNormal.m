function [estimate, samples, evaluations] = ISmonteCarloNormal(runNumber,numSteps,seed,mu1,sigma1,mu2,sigma2,distLimit,sumType)

fstepsX = generateRandomWalkStepsNormal(runNumber,numSteps,seed(1));
fstepsY = generateRandomWalkStepsNormal(runNumber,numSteps,seed(2));
fstepsZ = generateRandomWalkStepsNormal(runNumber,numSteps,seed(3));
frandomWalk = cat(3,fstepsX,fstepsY,fstepsZ);
fsamples = frandomWalk;

frandomWalkT = cat(3,(frandomWalk(:,:,1))',(frandomWalk(:,:,2))',(frandomWalk(:,:,3))');
fposEnd = sum(frandomWalkT);
% xyzEnd = squeeze(posEnd(:,end,:))
fposEnd = squeeze(fposEnd);
fdistEnd = vecnorm(fposEnd'); % If A is a matrix, then vecnorm returns the norm of each column.




% 
% fsteps = generateRandomWalkSteps(runNumber,numSteps,pThresholdf,seed);
% frandomWalk=generateRandomWalk(fsteps);
% %[fevaluation, fevaluationarray]=evaluatorFunc(frandomWalk,distLimitf);
% f_x = prob_x(fsteps,pThresholdf);


qstepsX = generateRandomWalkStepsNormalIS(runNumber,numSteps,mu2,sigma2,seed(1));
qstepsY = generateRandomWalkStepsNormalIS(runNumber,numSteps,mu2,sigma2,seed(2));
qstepsZ = generateRandomWalkStepsNormalIS(runNumber,numSteps,mu2,sigma2,seed(3));
qrandomWalk = cat(3,qstepsX,qstepsY,qstepsZ);
qsamples = qrandomWalk;

qrandomWalkT = cat(3,(qrandomWalk(:,:,1))',(qrandomWalk(:,:,2))',(qrandomWalk(:,:,3))');
qposEnd = sum(qrandomWalkT);
% xyzEnd = squeeze(posEnd(:,end,:))
qposEnd = squeeze(qposEnd);
qdistEnd = vecnorm(qposEnd'); % If A is a matrix, then vecnorm returns the norm of each column.

q_x=normalDistPDF(mu2,sigma2,qsamples);
f_x=normalDistPDF(mu1,sigma1,qsamples);


% qsteps = generateRandomWalkSteps(runNumber,numSteps,pThresholdq,seed);
% qrandomWalk=generateRandomWalk(qsteps);
% % [qevaluation, qevaluationarray]=evaluatorFunc(proposalFunc,distLimitq);
% q_x = prob_x(qsteps,pThresholdq);
% % plot(qrandomWalk','LineWidth',1.5,'Color',[0, 0.5, 0.5, 0.2])
% % hist(qrandomWalk(:,end))

% g_x = indFunc(qdistEnd,distLimit);
% g_xThreshold = 1;
g_x  = evaluatorFunc(qdistEnd,distLimit);
plikelihoodfq = zeros(runNumber,1);
likelihoodfq   = f_x./q_x;

for i=1:runNumber
likelihoodStep = squeeze(likelihoodfq(i,:,:));
pEachStep = prod(likelihoodStep');
pEachRun = prod(pEachStep);
plikelihoodfq(i,1) = pEachRun;
end

% plikelihoodfq  = prod(likelihoodfq');
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