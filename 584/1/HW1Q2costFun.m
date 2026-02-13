function F = HW1Q2costFun(x,L1,L2,R1,R2)
    f1 = ((x(1)-L1(1)).^2 + (x(2)-L1(2)).^2 - R1^2);
    f2 = ((x(1)-L2(1)).^2 + (x(2)-L2(2)).^2 - R2^2);
    F  = f1.^2+f2.^2;
end