
function evaluatorFuncOut = evaluatorFunc(totDist, threshold)

% totDist = randomWalks(:,end);
g_xTotDist = totDist;
g_xTotDist (totDist <= threshold) = 0;
g_xTotDist (totDist > threshold) = 1;
probDistLimit = sum(g_xTotDist)/length(g_xTotDist);
evaluatorFuncOut = g_xTotDist;
% evaluatorFuncOut = probDistLimit;
end