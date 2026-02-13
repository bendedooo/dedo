clear all; clc; close all;
%% part a
% a = 4/(4+pi);
% t = 4
Y = [];
% xx = @(x) a*sin(t);

fun = @(x,a) (a.*sin(x)).^2 + (a.*cos(x)).^2;
% L = (0.5)*integral(fun,0,pi/2);

% K = (xx-1)^2;

% 4/(4+pi);

J = @(a) (a*sin(pi/2) - 1).^2 + (0.5) * integral(@(x) fun(x,a),0,pi/2);
lim = 2;
avec = -lim:0.01:lim;
figure(1)
for i=1:length(avec)
a = avec(i);
y = J(a);
Y = [Y y];
end
hold on

% y = J(avec);
plot(avec,Y)
[Jmin,aind] = min(Y);
a_min = avec(aind)
J_min = Jmin
scatter(a_min,J_min)

a_ = 4/(4+pi)
J_a_ = J(a_)
%% part b
clear all; clc; close all;

syms a b x
K =  (a*pi/2 + b/2*(pi/2)^2 - 1)^2;
L =  (a*x + 0.5*b*x.^2).^2 + (a+b*x).^2;
J =  K+  0.5* int(L,[0 pi/2]) %bunlari fonksiyon yapmak istiyorum

dJ_da = diff(J,a)
dJ_db = diff(J,b)
[solx,soly] = solve(dJ_da == 0, dJ_db == 0);
a_min = double(solx)
b_min = double(soly)

clear x
% K = @(a,b) (a*pi/2 + b/2*(pi/2)^2 - 1)^2;
% L = @(a,b,x) (a*x + 0.5*b*x.^2).^2 + (a+b*x).^2;
K =  @(a,b) (a*pi/2 + b/2*(pi/2)^2 - 1)^2;
L =  @(a,b,x) (a*x + 0.5*b*x.^2).^2 + (a+b*x).^2;
J_min = K(a_min,b_min) + 0.5* integral(@(x) L(a_min,b_min,x),0, pi/2)
% J(a_min,b_min)


dJ_da = diff(J,a)
dJ_db = diff(J,b)

% check for strich convexity with the hessian matrix
dJ_da_da = diff(dJ_da,a);
dJ_da_db = diff(dJ_db,a);
dJ_db_da = diff(dJ_da,b);
dJ_db_db = diff(dJ_db,b);



HM = [dJ_da_da   dJ_da_db;...
      dJ_db_da   dJ_db_db]
eigHM = eig(HM)
isposdef = all(eigHM > 0)


% 
% % olmadi burasi hessian fonksiyonunu kullanarak a ve b min yerine koyarak yapmak istedim H olusturduktan
% % sonra ama a ve b sembolik oldugu icin substitute edemedim
% H = hessian(J,[a,b]);
% Hh = [matlabFunction(H(1,1)   matlabFunction(H(1,2);...
%       matlabFunction(H(2,1)   matlabFunction(H(2,2)];
% Hh(a_min,b_min);
% chol(Hh) %doesnt fail, pos def
% d = eig(Hh)
% isposdef = all(d > 0)
% double(d)

% bu error veriyo variable kalmadigi icin
% % a = a_min;
% % b = b_min;
% dJ_da_da = matlabFunction(dJ_da_da);
% dJ_da_db = matlabFunction(dJ_da_db);
% dJ_db_da = matlabFunction(dJ_db_da);
% dJ_db_db = matlabFunction(dJ_db_db);
% 
% HM = [dJ_da_da(a_min,b_min)   dJ_da_db(a_min,b_min);...
%       dJ_db_da(a_min,b_min)   dJ_db_db(a_min,b_min)]


% % fonskiyona cevirmeyi denedim
% hessianCheckConvexity(J,a_min,b_min)
% 
% function isConvex = hessianCheckConvexity(J,a_min,b_min)
% 
% J = matlabFunction(J)
% 
% dJ_da = diff(J,a)
% dJ_db = diff(J,b)
% 
% dJ_da_da = diff(dJ_da,a);
% dJ_da_db = diff(dJ_db,a);
% dJ_db_da = diff(dJ_da,b);
% dJ_db_db = diff(dJ_db,b);
% 
% % % we did not need to use the a* and b* here since second derivatives don't
% % % have variables
% % a = a_min;
% % b = b_min;
% 
% 
% HM = [dJ_da_da   dJ_da_db;...
%       dJ_db_da   dJ_db_db]
% eigHM = eig(HM)
% isposdef = all(eigHM > 0)
% % 
% % 
% % H = hessian(J,[a,b])
% % chol(H)
% % d = eig(H)
% % isposdef = all(d > 0)
% % double(d)
% 
% end