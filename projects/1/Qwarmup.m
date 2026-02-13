% Warmup deneme

clear all
close all
clc

seed = 1333; % rng seed for reproducibility


sumType = 0;
alpha_p = 3/2;
sampleNum = 1000;
n = 100;

MCest = zeros(1,n);

% for i=1:n
% % inverse sampling
% u = rand(1,1000);
% x = (1-u).^(-1/alpha_p);
% 
% % hist(x);
% MCEsti=sum(x)/sampleNum;
% MCest(1,i) = MCEsti;
% 
% end

sumType = 1;
[estimate, samples, evaluations] = monteCarlow(n,sampleNum,seed,alpha_p,sumType);


plot(estimate')
title('Warmup Question', 'Running Estimate')
xlabel('Sample Number')
ylabel('S_i')
legend('X_n','Location','best')
grid on

j=1;
for i=1:sampleNum
   
varn(j)=var(estimate(:,i));
 j=j+1;
end

semilogy(varn)

semilogy(varn)
hold on
semilogy(1./(1:sampleNum))
title('Warmup Question', 'Variance of Estimates')
xlabel('Sample Number')
ylabel('Variance')
legend('X_n','1/n','Location','best')
grid on
% runNumber = 100;
% [estimate, samples, evaluations] = monteCarlo1(runNumber,seed,alpha_p,sumType);

 
% if x >= 1
%     P_x = x^((alpha)/(alpha+1));
% elseif x<1
%     P_x = 0;
% end

% uniformDist = rand(1,1000)
% paretoDist = 1 ./ (uniformDist.^(1/alpha_p))
% hist(uniformDist)
% hist(paretoDist)


% syms x
% PDF = alpha_p/(x^(alpha_p+1));
% CDF = int(PDF);
% CDF_num= matlabFunction(CDF)


figure(3)
i=200:400:sampleNum;
h(1)=histogram(evaluations(:,i(1)),'FaceColor','b','FaceAlpha',0.1)
hold on
h(2)=histogram(evaluations(:,i(2)),'FaceColor','g','FaceAlpha',0.3)
h(3)=histogram(evaluations(:,i(3)),'FaceColor','r','FaceAlpha',0.5)
title('Warmup Question', 'Histogram of Evaluations')
xlabel('Sample Number')
ylabel('Frequency')
legend('n=200','n=600','n=1000','Location','best')
grid on

figure(4)
i=200:400:sampleNum;
h(1)=histogram(estimate(:,i(1)),'FaceColor','b','FaceAlpha',0.1)
hold on
h(2)=histogram(estimate(:,i(2)),'FaceColor','g','FaceAlpha',0.3)
h(3)=histogram(estimate(:,i(3)),'FaceColor','r','FaceAlpha',0.5)
title('Warmup Question', 'Histogram of Estimates')
xlabel('Sample Number')
ylabel('Frequency')
legend('n=200','n=600','n=1000','Location','best')
grid on

z = 1.96;
stdErrorrunning = std(estimate)./sqrt(n);
var_ = var(estimate(:,end));
intervalrunning = [mean(estimate(:,end))-z*stdErrorrunning(end); mean(estimate(:,end))+z*stdErrorrunning(end)];
confP= sum((intervalrunning(1) <= estimate(:,end)) & (estimate(:,end) <= intervalrunning(2)))/n;

