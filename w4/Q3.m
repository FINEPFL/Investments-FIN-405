%% Homework 4 - Investments
% Exercise 3 - Evaluating the momentum strategy
clear all;
close all;
clc;

Alldata = xlsread('MOM');
% by order: MKT SMB HML MOM
%% (a)
% For period 1 going from January 1927 to latest date
AllData = Alldata(:, 2:5); 
Mean1 = mean(AllData);
Volatility1 = std(AllData);
SharpeRatio1 = Mean1./Volatility1;
[h1,p1,ci1,stats1] = ttest(AllData);

% For period 2 going from July 1963 to latest date (439:end,:)

Mean2 = mean(AllData(439:end,:));
Volatility2 = std(AllData(439:end,:));
SharpeRatio2 = Mean2./Volatility2;
[h2,p2,ci2,stats2] = ttest(AllData(439:end,:));

% For period 3 going from January 1991 to latest date (769:end,:)

Mean3 = mean(AllData(769:end,:));
Volatility3 = std(AllData(769:end,:));
SharpeRatio3 = Mean3./Volatility3;
[h3,p3,ci3,stats3] = ttest(AllData(769:end,:));


%% (b)
% For period 1
Positive1 = sign(AllData);
Positive1( Positive1==-1 )=0;
Fraction1 = sum(Positive1)./length(Positive1);

% For period 2
Positive2 = sign(AllData(439:end,:));
Positive2( Positive1==-1 )=0;
Fraction2 = sum(Positive2)./length(Positive2);

% For period 3
Positive3 = sign(AllData(769:end,:));
Positive3( Positive1==-1 )=0;
Fraction3 = sum(Positive3)./length(Positive3);


%% (c) Portfolio combination
ExcessRetunrNewStrategy = 0.5*(AllData(:,3)+AllData(:,4)); % Denotes as strategy 5
MeanNewStrategy = mean(ExcessRetunrNewStrategy);
VolatilityNewStrategy = std(ExcessRetunrNewStrategy);
SharpeRatioNewStrategy = MeanNewStrategy./VolatilityNewStrategy;
[hN,pN,ciN,statsN] = ttest(ExcessRetunrNewStrategy);

%% (d)
AllData(:,5) = ExcessRetunrNewStrategy;
W = [1,1,1,1,1];
%W(1,:) = 1;
for j=1:size(AllData,2)
    for i=1:length(AllData)
        W(i+1,j) = W(i,j)*(1+AllData(i,j));
    end
end

figure
plot(Date, W(2:end, 1), '.', Date, W(2:end, 2), '-', Date, W(2:end, 3), '--', Date, W(2:end, 4), '-.', 'linewidth', 1.5);
xlabel('Date')
ylabel('$ amount')
legend('Market','SMB','HML','MOM','50% HML + 50% UMD','Interpreter','Latex','Location','Best', 'fontsize', 20)
