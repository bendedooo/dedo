function psuedoInvOut = pseudoInvFunc(AMatrix)

psuedoInvOut= inv(AMatrix'*AMatrix)*AMatrix';

end