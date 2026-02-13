function F = HW2Q2costFun(x,L1,L2,theta,psi,XS)
% f1 = (x(1)-L1(1)) * (L2(1)-L1(1)) ...
%     + (x(2)-L1(2)) * (L2(2)-L1(2)) ...
%     + sqrt( (x(1)-L1(1))^2 + (x(2)-L1(2))^2 ) * sqrt( (x(1)-L2(1))^2 + (x(2)-L2(2))^2 ) * cos(theta)...
%     - ((x(1)-L1(1))^2+(x(2)-L1(2))^2);
% f2 = ( (L2(1)-x(1))*0 + (L2(2)-x(2))*1 ) ...
%     - sqrt((x(1)-L2(1))^2 + (x(2)-L2(2))^2) * cos(psi);

PL1  = x-L1;
PL1  = PL1';
L2P  = L2-x;
L2P  = L2P';
L2L1 = L2-L1;
L2L1 = L2L1';
% XS   = [0 1 0]';

f1 = PL1'* L2L1 + norm(PL1)*norm(-L2P)*cos(theta)-norm(PL1)^2;
f2 = L2P'* XS - norm(-L2P)*cos(psi);

F  = f1^2+f2^2;
end