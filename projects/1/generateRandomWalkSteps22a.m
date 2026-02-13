% generate a set of steps for a random walk

function X = generateRandomWalkStepsNormal(runNumber,numSteps)
% rng seed fix
% s = rng;
% s.Seed = seed;
% rng(s.Seed)
% for i=1:numSteps
% xi = randn(1,3);
% norm(xi)
% end
X = randn(runNumber,numSteps); % sample from a unifrom
end