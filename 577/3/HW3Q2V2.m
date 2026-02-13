%ncdisp('sst.wkmean.1990-present.nc')
clear; close all; clc
% datpath = '../DATA/';
figpath = '/Users/dedo/Documents/UMICH/Winter 22/AEROSP577/HW/3/kod/';
% seed random number generator for reproducibility
% rng(729);
s = rng;
s.Seed = 729;
rng(s.Seed);
%% Part a 
[ Lat, Lon, time, mask, sst ] = read_data_enso( 'sst.wkmean.1990-present.nc','lsmask.nc');

% each element of time array is a new week, in units of days
t0 = datetime(1800,1,1,0,0,0) + days(time(1));
tfin = datetime(1800,1,1,0,0,0) + days(time(end));

[m,n,p] = size(sst);
N = m*n;
X = zeros(N,p);
M = length(mask(mask==1));
Y = zeros(M,length(time));
for i=1:length(time)
   snapshot = reshape(sst(:,:,i),N,1);   
   Y(:,i) = snapshot(mask==1);
end

% train on first 600 snapshots
Itrain = 1:600; 
Itest = 800; 

TrainRaw = Y(:,Itrain);
timeavg = mean(TrainRaw,2);    % M = mean(A,dim) returns the mean along dimension dim. For example, if A is a matrix, 
                               % then mean(A,2) is a column vector containing the mean of each row.
Train = bsxfun(@minus,TrainRaw,timeavg);  % Train = Train - Time Average

[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);
thresh = optimal_SVHT_coef(n/m,0)*median(sing);
r_opt = length(sing(sing>=thresh));

R = [100 200 r_opt 600];

figure(1)
plot(sing,'.','color',.7*[1 1 1]);
hold on
plot(sing(1:r_opt),'b.');

plot(R, sing(R),'ro');
plot(R, sing(R),'r.');
set(gca,'yscale','log','xlim',[0 length(sing)+1]);
grid on

set(gca,'ylim',[1 10^5]);
ylabel('\sigma_{snapshot}')
xlabel('Snapshots')
title('Singular Value Spectrum')

cumEnergy = cumsum(sing)./sum(sing); % normalized cumulative energy = normalized sum of singular values 
figure(2)
plot(cumEnergy, LineWidth=2);
grid on
ylabel('Cumulative Energy')
xlabel('Snapshots')
title('Cumulative Energy')
%% Part b
close all
meansst = timeavg;
% select validation snapshot
x = Y(:,Itest(1))-timeavg;
bounds = [min(x+timeavg) max(x+timeavg)];


r = 150; 
% POD approximation with r eigenssts
PsiNew = Psi(:,1:r);
xproj = PsiNew*(PsiNew'*x);

figure(2)
display_fig(xproj+timeavg,mask,[],bounds);
title('800^{th} Snapshot Projected on the New Basis')
figure(3)
display_fig(x+timeavg,mask,[],bounds);
title('800^{th} True Snapshot')
%% Part c
 
% QDEIM with r QR sensors
r = 50;
[~,~,pivot] = qr(PsiNew','vector');
sensors = pivot(1:r);

figure(4)
display_sensors(x+meansst,mask,sensors);
title('QR Sampling Locations')

xls = PsiNew*(PsiNew(sensors,:)\x(sensors));
figure(5)
display_fig(xls+meansst,mask,[],bounds);
title('800^{th} Snapshot Reconstructed with Pivoted QR Decompoition (50 Sensors)')

errorQR = norm(xls-x)/norm(x)
errorQR = norm(xls-x)

%% Part d
% QDEIM with r QR sensors
r = 150;
% [~,~,pivot2] = qr(Psi(:,1:r)','vector');
sensors = pivot(1:r);

figure(6)
display_sensors(x+meansst,mask,sensors);
title('QR Sampling Locations')

xls = PsiNew*(PsiNew(sensors,:)\x(sensors));
figure(7)
display_fig(xls+meansst,mask,[],bounds);
title('800^{th} Snapshot Reconstructed with Pivoted QR Decompoition (150 Sensors)')

errorQR = norm(xls-x)/norm(x)
errorQR = norm(xls-x)

%% Part e
% QDEIM with r QR sensors
% sensind = randperm(m,ceil(length(PsiNew)/4));
sensind = 1:4:m;
PsiNew2=PsiNew(sensind,:);
% PsiNew2=imresize(PsiNew,[ceil(length(PsiNew)/4),150]);
r = 250;
[Q,R,pivot] = qr(PsiNew2*PsiNew2','vector');
pivot = sensind(pivot);
sensors = pivot(1:r);

figure(8)
display_sensors(x+meansst,mask,sensors);
title('QR Sampling Locations')

xls = PsiNew*(PsiNew(sensors,:)\x(sensors));
figure(9)
display_fig(xls+meansst,mask,[],bounds);
title('800^{th} Snapshot Reconstructed with Pivoted QR Decompoition (250 Sensors)')

errorQR = norm(xls-x)/norm(x)
norm(xls-x)
%% Part f
% Random reconstruction with r sensors
R = [50 150 250];

for ii = 1:length(R)
    r = R(ii);
    sensors = randperm(m,r);

    figure((ii+4)*2)
    display_sensors(x+meansst,mask,sensors);
    title('Random Sampling Locations')

    figure((ii+4)*2+1)
    xls = PsiNew*(pinv(PsiNew(sensors,:))*x(sensors));
    display_fig(xls+meansst,mask,[],bounds);
    title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L2 Min',r))

    errorRandomL2 = norm(xls-x)/norm(x)
    norm(xls-x)
end
%% Part f
R = [50 150 250];

    % first obtain sensors
sens1 = randperm(m,R(1));
sens2 = randperm(m,R(2));
sens3 = randperm(m,R(3));

% form large Theta matrices all at once
Theta1 = zeros(R(1),n);
Theta2 = zeros(R(2),n);
Theta3 = zeros(R(3),n);
Isp = speye(m);
% Isp = eye(m);

for i = 1:m
    psi = idct2(Isp(:,i));
    Theta1(:,i) = psi(sens1);
    Theta2(:,i) = psi(sens2);
    Theta3(:,i) = psi(sens3);

end
% r = 50 sensors
r = R(1);
y = x(sens1);%+meansst(sens1);

cvx_begin quiet;
variable a(m);
minimize(norm(a,2)) ;
subject to 
Theta1*a  == y;
cvx_end;

xcs = idct2(a);
figure(16) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L2 Min',r))
errorRandomL150 = norm(xcs-x)/norm(x)
norm(xcs-x)

% r = 150 sensors
r = R(2);
y = x(sens2);%+meansst(sens2);

cvx_begin quiet;
variable a(m);
minimize(norm(a,2)) ;
subject to 
Theta2*a  == y;
cvx_end;

xcs = idct2(a);
figure(17) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L2 Min',r))
errorRandomL1150 = norm(xcs-x)/norm(x)
norm(xcs-x)


% r = 250 sensors
r = R(3);
y = x(sens3);%+meansst(sens3);

cvx_begin quiet;
variable a(m);
minimize(norm(a,2)) ;
subject to 
Theta3*a  == y;
cvx_end;

xcs = idct2(a);
figure(18) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L2 Min',r))
errorRandomL1250 = norm(xcs-x)/norm(x)
norm(xcs-x)

%% Part g
R = [50 150 250];

    % first obtain sensors
sens1 = randperm(m,R(1));
sens2 = randperm(m,R(2));
sens3 = randperm(m,R(3));

% form large Theta matrices all at once
Theta1 = zeros(R(1),n);
Theta2 = zeros(R(2),n);
Theta3 = zeros(R(3),n);
% Isp = speye(m);
Isp = eye(m);

for i = 1:m
    psi = idct(Isp(:,i));
    Theta1(:,i) = psi(sens1);
    Theta2(:,i) = psi(sens2);
    Theta3(:,i) = psi(sens3);

end
% r = 50 sensors
r = R(1);
y = x(sens1);%+meansst(sens1);

cvx_begin quiet;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta1*a  == y;
cvx_end;

xcs = idct(a);
figure(16) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L1 Min',r))
errorRandomL150 = norm(xcs-x)/norm(x)
norm(xcs-x)
norm((xcs-x),1)/norm(x,1)

% r = 100 sensors
r = R(2);
y = x(sens2);%+meansst(sens2);

cvx_begin quiet;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta2*a  == y;
cvx_end;

xcs = idct(a);
figure(17) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L1 Min',r))
errorRandomL1100 = norm(xcs-x)/norm(x)
norm(xcs-x)
norm((xcs-x),1)/norm(x,1)

% r = 150 sensors
r = R(3);
y = x(sens3);%+meansst(sens3);

cvx_begin quiet;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta3*a  == y;
cvx_end;

xcs = idct(a);
figure(18) 
display_fig(xcs+meansst,mask,[],bounds);
title(sprintf('800^{th} Snapshot Reconstructed with {%d} Random Sensors L1 Min',r))
errorRandomL1150 = norm(xcs-x)/norm(x)
norm(xcs-x)
norm((xcs-x),1)/norm(x,1)