% q =  [sin(0.6)*[1 0 0] cos(0.6)]
% qb = [q(4); q(1); q(2); q(3)]

for ii = 1:length(timeVec)
    q = Xorb(ii,10:13)';
    a = quaToDCM(q);
    DCMmult = a*a';
    b = quaToDCM(q);
    DCMmultb = b*b';
    % b = qToDCM(q)
    isOrth(ii,1:9) = reshape(DCMmult,[1 9]);
    isOrthb(ii,1:9) = reshape(DCMmultb,[1 9]);

    isqCorr(ii)=q'*q;
end