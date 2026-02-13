
function probWalk = prob_xNormal(inputArray,mu,sigma)

probWalk = normalDistPDF(mu,sigma,inputAray);
% 
% inputArray (inputArray ==  -1) = pValue ;
% inputArray (inputArray ==  1) = 1-pValue ;
% probWalk = (inputArray);
end