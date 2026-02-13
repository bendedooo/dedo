
function probWalk = prob_x(inputArray,pValue)
inputArray (inputArray ==  -1) = pValue ;
inputArray (inputArray ==  1) = 1-pValue ;
probWalk = (inputArray);
end