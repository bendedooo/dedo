function [estimate, samples, evaluations] = monteCarlo2(runNumber,numSteps,seed,pThreshold,distLimit,sumType)

steps = generateRandomWalkSteps(runNumber,numSteps,pThreshold,seed);
randomWalk=generateRandomWalk(steps);
samples = randomWalk;
totDist = randomWalk(:,end);
evaluations =evaluatorFunc(totDist,distLimit);



% samples = sampleGenerator(numSamples);
% evaluations = gEvaluate(samples);
% sumType
if sumType == 0
    estimate = sum(evaluations) / (runNumber); % indicator func sum columnwise summation
else
    estimate = cumsum(evaluations) / (1:runNumber); % columnwise summation
end
end