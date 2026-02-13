function [estimate, samples, evaluations] = monteCarloNormalIS(runNumber,numSteps,seed,distLimit,sumType)

stepsX = generateRandomWalkStepsNormalIS(runNumber,numSteps,seed(1));
stepsY = generateRandomWalkStepsNormalIS(runNumber,numSteps,seed(2));
stepsZ = generateRandomWalkStepsNormalIS(runNumber,numSteps,seed(3));
randomWalk = cat(3,stepsX,stepsY,stepsZ);
samples = randomWalk;

randomWalkT = cat(3,(randomWalk(:,:,1))',(randomWalk(:,:,2))',(randomWalk(:,:,3))');
posEnd = sum(randomWalkT);
% xyzEnd = squeeze(posEnd(:,end,:))
posEnd = squeeze(posEnd);
distEnd = vecnorm(posEnd'); % If A is a matrix, then vecnorm returns the norm of each column.


% posEnd = samples(:,end,:);
% xyzEnd = squeeze(posEnd);
% distEnd = vecnorm(xyzEnd'); % If A is a matrix, then vecnorm returns the norm of each column.

[evaluations, gtot]=evaluatorFunc(distEnd,distLimit);


if sumType == 0
    estimate = sum(evaluations) / (runNumber); % indicator func sum columnwise summation
else
    estimate = cumsum(evaluations) / (1:runNumber+1); % columnwise summation
end
end