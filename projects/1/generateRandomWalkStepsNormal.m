% generate a set of steps for a random walk

function X = generateRandomWalkStepsNormal(runNumber,numSteps,seed)
% rng seed fix
s = rng;
s.Seed = seed;
rng(s.Seed)

X = randn(runNumber,numSteps); % sample from a unifrom

end