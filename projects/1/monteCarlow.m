function [estimate, samples, evaluations] = monteCarlow(runNumber,numSteps,seed,alpha_p,sumType)
s = rng;
s.Seed = seed;
rng(s.Seed);

u = rand(runNumber,numSteps);
u (u>0.9) = 0.1;
x = (1-u).^(-1/alpha_p);
samples = x;
evaluations = x;

% samples = sampleGenerator(numSamples);
% evaluations = gEvaluate(samples);
% sumType
if sumType == 0
    estimate = sum(evaluations') / (numSteps); % indicator func sum columnwise summation
else
    indMatrix = ones(runNumber,1)*(1:numSteps);
    cumsumMatrix = cumsum(evaluations');
    estimate = cumsumMatrix' ./ (1:numSteps); % columnwise summation
end
end