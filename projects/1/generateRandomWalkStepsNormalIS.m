% generate a set of steps for a random walk

function X = generateRandomWalkStepsNormalIS(runNumber,numSteps,mu,sigma,seed)
% rng seed fix
s = rng;
s.Seed = seed;
rng(s.Seed)

X = sigma*randn(runNumber,numSteps) + mu; % sample from a unifrom

end