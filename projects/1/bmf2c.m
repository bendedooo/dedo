function bmf2cOut = bmf2c(tFinal, dt, runNumber, bf, m)
sampleNum = ceil(tFinal/dt)+1;
bc = zeros(sampleNum,runNumber);
bc(1,:) = 0;
for i = 2:sampleNum
    delta = bf((i-1)*m+1,:) - bf((i-2)*m+1,:);
    bc(i,:) = bc(i-1,:) + delta;
end
bmf2cOut = bc;
end