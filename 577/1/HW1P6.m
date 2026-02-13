clc
clear all
close all

%% import data and divide into equal batches to create new datasets
givenData = importdata('hw1data.txt');
tData = givenData(:,1);
yData = givenData(:,2);

indexNum = size(givenData,1);
% IndexOrder = 1:indexNum;
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

% tData1 = [tDataBatches(:,2); tDataBatches(:,3)];
% yData1 = [yDataBatches(:,2); yDataBatches(:,3)];
%
% tData2 = [tDataBatches(:,1); tDataBatches(:,3)];
% yData2 = [yDataBatches(:,1); yDataBatches(:,3)];

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

    % r = y - A(x_hat) bunu mu istiyo validation error ne???
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


d = [2 3 4 5];
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

%% plot
faceColorMatrix = [0 .7 .7;
    0.1 0.7 0.1;
    0.7 0.7 0.1;
    0.3 0.1 0.1];

figure(1)
t = sort(tDataBatches(:,4));
tVec = t(1):0.01:t(end);
tVec = tVec';

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
%% cop
% %% plot
% faceColorMatrix = [0 .7 .7;
%     0.1 0.7 0.1;
%     0.7 0.7 0.1;
%     0.3 0.1 0.1];
%
% figure(1)
% for ii = 1:length(d)
%     dii = d(ii);
%     MinMSEInd = find(MSE(:,1)==min(MSE(:,1))); % d=ii (2,3,4,5)
%     xii = X(1:dii,MinMSEInd,ii);
%     Aii = [];
%     for jj= 1:dii
%         Aii = [Aii sort(tDataBatches(:,4)).^(jj-1)];
%     end
%     yii = Aii*xii;
%     plot(sort(tDataBatches(:,4)), yii,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
%     hold on
% end
%
% scatter(tDataBatches(:,4), yDataBatches(:,4),'MarkerEdgeColor',[.5 0.7 .9],'MarkerFaceColor',[.5 0.7 .9],'LineWidth',1)
%
% xlabel('t')
% ylabel('y')
%
%
% legend('d=2','d=3','d=4','d=5','Data Batch 4')
% title('Problem 4: Part A.4 Best Model Predictions with "d" Values')

% d2MinMSEInd = find(MSE(:,1)==min(MSE(:,1))); % d=2
% d3MinMSEInd = find(MSE(:,2)==min(MSE(:,2))); % d=3
% d4MinMSEInd = find(MSE(:,3)==min(MSE(:,3))); % d=4
% d5MinMSEInd = find(MSE(:,4)==min(MSE(:,4))); % d=5
%% part b

lambda = 0.33; %????

d = 5;
A =[];

for jj= 1:d
    A = [A (tData).^(jj-1)];
end

b = yData;
n = d; % ????

cvx_begin quiet
variable x(n)
minimize( square_pos(norm(A*x-b,2))+ lambda*square_pos(norm(x,2)) )
cvx_end
cvx_status

figure(2)
% t = sort(tData);
% tVec = t(1):0.01:t(end);
% tVec = tVec';
tVec = tData;
dii=5;

Aii = [];

for jj= 1:dii
    Aii = [Aii tVec.^(jj-1)];
end
yii = Aii*x;
plot(tVec,yii,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
hold on

tDataPlot = sort(tData);
scatter(tData, yData,'MarkerEdgeColor',[.5 0.7 .9],'MarkerFaceColor',[.5 0.7 .9],'LineWidth',1)

xlabel('t')
ylabel('y')
legend('d=5','Data')
title('Problem 6: Part B')

xCVX = x;
xAnalytical = ( inv(Aii'*Aii + lambda*eye(dii)) * Aii')*yData;
%% part c
lambdaArray = 0:0.05:2;
d = [2 3 4 5];

newIndexOrder = randperm(indexNum);

tDataTrain = tData(newIndexOrder(1:indexNum*0.8));
tDataValid = tData(newIndexOrder(indexNum*0.8+1:end));

yDataTrain = yData(newIndexOrder(1:indexNum*0.8));
yDataValid = yData(newIndexOrder(indexNum*0.8+1:end));

% tDataTrain = tData(1:indexNum*0.8);
% tDataValid = tData(indexNum*0.8+1:end);
% 
% yDataTrain = yData(1:indexNum*0.8);
% yDataValid = yData(indexNum*0.8+1:end);

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
        minimize( (square_pos(norm(ATrain*x-b,2)))+ (lambda*square_pos(norm(x,2))) )
        cvx_end
        cvx_status

        yValid = AValid*x;
        rValid = yValid-yDataValid;
        MSE = sum(rValid.^2)/length(rValid);
        MSEArray = [MSEArray MSE];
    end
        figure(3)

        plot(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
end

xlabel('\lambda')
ylabel('MSE')


legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part C')

%% part d
% lambdaArray = 0:0.05:2;
% d = [2 3 4 5];
% 
% tDataTrain = tData(1:indexNum*0.8);
% tDataValid = tData(indexNum*0.8+1:end);
% 
% yDataTrain = yData(1:indexNum*0.8);
% yDataValid = yData(indexNum*0.8+1:end);

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
        cvx_status

        yValid = AValid*x;
        rValid = yValid-yDataValid;
        MSE = sum(rValid.^2)/length(rValid);
        MSEArray = [MSEArray MSE];
    end
        figure(4)

        plot(lambdaArray,MSEArray,'LineWidth',1.5,'Color',faceColorMatrix(ii,:));
        hold on
end

xlabel('\lambda')
ylabel('MSE')


legend('d=2','d=3','d=4','d=5')
title('Problem 6: Part D')
