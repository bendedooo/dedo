function out = CMatrix(posEst,posL)
% posL = colums ith L
g = sqrt((posEst(1,1)-posL(1,:)).^2+(posEst(2,1)-posL(2,:)).^2);
out = [(posEst(1,1)-posL(1,:))' (posEst(2,1)-posL(2,:))']./g';
end