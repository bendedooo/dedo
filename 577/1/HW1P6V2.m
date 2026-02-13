clc
clear all
close all
%% import data and divide into equal batches to create new datasets
givenData = importdata('hw1data.txt');
tData = givenData(:,1);
yData = givenData(:,2);

indexNum = size(givenData,1);
newIndexOrder = randperm(indexNum);
batchNum = 4;
batchSize = floor(indexNum/batchNum);
tDataBatches = zeros(batchSize,batchNum);
yDataBatches = zeros(batchSize,batchNum);
batchIndexMatrix = zeros(batchSize,batchNum);

for i = 1:batchNum %floor(indexNum/batchSize);
    batchIndexMatrix(:,i) = newIndexOrder(batchSize*(i-1)+1:batchSize*i);
    tDataBatches(:,i) = tData(batchIndexMatrix(:,i));
    yDataBatches(:,i) = yData(batchIndexMatrix(:,i));
end

tDataNew = [];
yDataNew = [];
for i=1:3
    for j=1:3
        if i<j || i==~j
            tDataNew = cat(2,tDataNew,[tDataBatches(:,i); tDataBatches(:,j)]);
            yDataNew = cat(2,yDataNew,[yDataBatches(:,i); yDataBatches(:,j)]);
        end
    end
end
%% solve
d = [2 3 4 5];
MSE3 = [];
X = [];
for ii = 1:length(d)
    dii = d(ii);
    Aii1 = [];
    Aii2 = [];
    Aii3 = [];
    for jj= 1:dii
        Aii1 = [Aii1 (tDataNew(:,1)).^(jj-1)];
        Aii2 = [Aii2 (tDataNew(:,2)).^(jj-1)];
        Aii3 = [Aii3 (tDataNew(:,3)).^(jj-1)];
    end
    % r = y - A(x_hat) 
    Aii1PI = pseudoInvFunc(Aii1);
    xii1 = Aii1PI*(yDataNew(:,1));
    rii1 = (yDataNew(:,1))-Aii1*Aii1PI*(yDataNew(:,1));
    MSEii1 = sum(rii1.^2)/length(rii1);
    Aii2PI = pseudoInvFunc(Aii2);
    xii2 = Aii2PI*(yDataNew(:,2));
    rii2 = (yDataNew(:,2))-Aii2*Aii2PI*(yDataNew(:,2));
    MSEii2 = sum(rii2.^2)/length(rii2);
    Aii3PI = pseudoInvFunc(Aii3);
    xii3 = Aii3PI*(yDataNew(:,3));
    rii3 = (yDataNew(:,3))-Aii3*Aii3PI*(yDataNew(:,3));
    MSEii3 = sum(rii3.^2)/length(rii3);
    xii = [xii1 xii2 xii3; zeros(d(end)-d(ii),1) zeros(d(end)-d(ii),1) zeros(d(end)-d(ii),1)];
    X   = cat(3,X,xii);
    MSEii = [MSEii1; MSEii2; MSEii3];
    MSE3 = [MSE3 MSEii];
end
%% part 4
MSE4 = [];
Y = [];
for ii = 1:length(d)
    dii = d(ii);
    Aii = [];
    for jj = 1:dii
        Aii = [Aii tDataBatches(:,4).^(jj-1)];
    end
    yii1 = Aii*X(1:dii,1,ii);
    rii1 = tDataBatches(:,4)-yii1;
    MSEii1 = sum(rii1.^2)/length(rii1);
    yii2 = Aii*X(1:dii,2,ii);
    rii2 = tDataBatches(:,4)-yii2;
    MSEii2 = sum(rii2.^2)/length(rii2);
    yii3 = Aii*X(1:dii,3,ii);
    rii3 = tDataBatches(:,4)-yii3;
    MSEii3 = sum(rii3.^2)/length(rii3);
    yii = [yii1 yii2 yii3];
    Y   = cat(3,Y,yii);
    MSEii = [MSEii1; MSEii2; MSEii3];
    MSE4 = [MSE4 MSEii];
end
faceColorMatrix = [0 .7 .7;0.1 0.7 0.1;0.7 0.7 0.1;0.3 0.1 0.1];
figure(1)
t = sort(tDataBatches(:,4));
tVec = t(1):0.01:t(end);
tVec = tVec';

% tVec = sort(tData);

for ii = 1:length(d)
    dii = d(ii);
    MinMSEInd = find(MSE4(:,1)==min(MSE4(:,1))); % d=ii (2,3,4,5)
    Aii = [];
    xii = X(1:dii,MinMSEInd,ii);
    for jj= 1:dii
        Aii = [Aii tVec.^(jj-1)];
    end
    yii = Aii*xii;
    plot(tVec,yii,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
    hold on
end
scatter(tDataBatches(:,4), yDataBatches(:,4),'MarkerEdgeColor',[.5 0.7 .9],'MarkerFaceColor',[.5 0.7 .9],'LineWidth',1)
xlabel('t')
ylabel('y')
legend('d=2','d=3','d=4','d=5','Data Batch 4')
title('Problem 6: Part A.4 Best Model Predictions with "d" Values')
%% part b
lambda = 0.33; 
dValue = 5;
A =[];
for jj= 1:dValue
    A = [A (tData).^(jj-1)];
end
b = yData;
n = dValue; 
cvx_begin quiet
variable x(n)
minimize( square_pos(norm(A*x-b,2))+ lambda*square_pos(norm(x,2)) )
cvx_end
cvx_status

figure(2)
tDataPlot = sort(tData);
yii = A*x;
plot(tDataPlot,yii,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
hold on
scatter(tData, yData,'MarkerEdgeColor',[.5 0.7 .9],'MarkerFaceColor',[.5 0.7 .9],'LineWidth',1)
xlabel('t')
ylabel('y')
legend('d=5','Data')
title('Problem 6: Part B')

xCVX = x;
xAnalytical = ( inv(A'*A + lambda*eye(dii)) * A')*yData;
%% part c
lambdaArray = linspace(10^-5,10,500);
% lambdaArray = logspace(-8,2,10);

d = [2 3 4 5];

newIndexOrder = randperm(indexNum);
tDataTrain = tData(newIndexOrder(1:indexNum*0.8));
tDataValid = tData(newIndexOrder(indexNum*0.8+1:end));
yDataTrain = yData(newIndexOrder(1:indexNum*0.8));
yDataValid = yData(newIndexOrder(indexNum*0.8+1:end));

for ii = 1:length(d)
    dii = d(ii);
    ATrain = [];
    AValid = [];
    for jj= 1:dii
        ATrain = [ATrain (tDataTrain).^(jj-1)];
    end
    for jj = 1:dii
%         AValid = [AValid tDataValid.^(jj-1)];
        AValid(:,jj)=tDataValid.^(jj-1);
    end
    MSEArray = [];
    for k=1:length(lambdaArray)
        lambda = lambdaArray(k);
        b = yDataTrain;
        n = dii; 
        cvx_begin quiet
        variable x(n)
        minimize( (0.5)*(square_pos(norm(ATrain*x-b,2)))+ (0.5)*(lambda*square_pos(norm(x,2))) )
        cvx_end
%         cvx_status
        yValid = AValid*x;
        rValid = yValid-yDataValid;
        MSE = sum(rValid.^2)/length(rValid);
        MSEArray = [MSEArray MSE];
    end
        figure(3)
        plot(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
        figure(4)
        semilogx(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
end
xlabel('log \lambda')
ylabel('MSE')
legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part C')
figure(3)
xlabel('\lambda')
ylabel('MSE')
legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part C')
%% part d
for ii = 1:length(d)
    dii = d(ii);
    ATrain = [];
    AValid = [];
    for jj= 1:dii
        ATrain = [ATrain (tDataTrain).^(jj-1)];
    end
    for jj = 1:dii
%         AValid = [AValid tDataValid.^(jj-1)];
        AValid(:,jj)=tDataValid.^(jj-1);
    end
    MSEArray = [];
    for k=1:length(lambdaArray)
        lambda = lambdaArray(k);
        b = yDataTrain;
        n = dii; % ????
        cvx_begin quiet
        variable x(n)
        minimize( square_pos((norm(ATrain*x-b,2)))+ lambda*(norm(x,1)) )
        cvx_end
%         cvx_status
        yValid = AValid*x;
        rValid = yValid-yDataValid;
        MSE = sum(rValid.^2)/length(rValid);
        MSEArray = [MSEArray MSE];
    end
        figure(6)
        plot(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
        figure(7)
        semilogx(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
end
xlabel('log \lambda')
ylabel('MSE')
legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part D')
figure(6)
xlabel('\lambda')
ylabel('MSE')
legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part D')