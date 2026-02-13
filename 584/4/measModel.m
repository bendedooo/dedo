function out = measModel(posEst,posL)
% posL = colums ith L
out = sqrt((posEst(1,1)-posL(1,:)).^2+(posEst(2,1)-posL(2,:)).^2)';
end