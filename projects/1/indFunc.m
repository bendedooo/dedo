
function indFuncOut = indFunc(inputArray,distLimit)
inputArray (inputArray <= distLimit) = 0;
inputArray (inputArray >  distLimit) = 1;
indFuncOut = inputArray;
end