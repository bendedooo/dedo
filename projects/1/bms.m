function bmsOut = bms(tFinal, dt, runNumber)


sampleNum = ceil(tFinal/dt);
samples = [zeros(runNumber,1) randn(runNumber,sampleNum)*sqrt(dt)];
bmsOut=cumsum(samples'); %If A is a matrix, then cumsum(A) returns a matrix containing the cumulative sums for each column of A.
end