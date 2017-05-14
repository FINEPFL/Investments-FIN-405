clearvars; close all; clc
raw_data = xlsread('Commodity_prices.xlsx');
% Crude Oil | Copper | Live Cattle | Cotton | Soybean | Live Hogs | 
% Sugar | Gold | Silver | Coffee | Wheat | Corn

% q1
A = raw_data(1:end-1, 2:2:24);
B = raw_data(2:end, 1:2:23);
returns = (B - A)./A;

ari_mean = mean(returns);
ret_std = std(returns);

cumed = cumprod((returns + 1), 1);
geo_mean = cumed(end, :).^(1/(size(cumed, 1)-1)) - 1;
skew = skewness(returns);
kurt = kurtosis(returns);
% geomean(returns+1)-1

% q2
corr_mat = corr(returns);
avg_corr = (sum(sum(corr_mat)) - sum(diag(corr_mat)))/2/66;

returns_1  = returns(1:0.5*size(returns, 1), :);
corr_mat_1 = corr(returns_1);
avg_corr_1 = (sum(sum(corr_mat_1)) - sum(diag(corr_mat_1)))/2/66;

returns_2 = returns(0.5*size(returns, 1)+1:end, :);
corr_mat_2 = corr(returns_2);
avg_corr_2 = (sum(sum(corr_mat_2)) - sum(diag(corr_mat_2)))/2/66;

close all

% q3
name_cell = {'crude oil', 'copper', 'live cattle', 'cotton', 'soybean', 'live hogs', 'suger' , 'gold', 'silver', 'coffee', 'wheat', 'corn'};
exce_return = mean(returns(0.5*size(returns, 1)+1:end, :));
roll_return = mean((-raw_data(0.5*size(raw_data, 1)+1:end, 2:2:24) + raw_data(0.5*size(raw_data, 1)+1:end, 1:2:23))./raw_data(0.5*size(raw_data, 1)+1:end, 2:2:24));
% exce_return = mean(returns);
% roll_return = mean((-raw_data(:, 2:2:24) + raw_data(:, 1:2:23))./raw_data(:, 2:2:24));

spot_return = exce_return - roll_return;
figure
for i=1:length(roll_return)
    plot(roll_return(i) * 12, exce_return(i) * 12, '.', 'color', [rand rand rand], 'markersize', 30); 
    hold on;
    text_hd = text(roll_return(i) * 12, exce_return(i) * 12, name_cell(i));
    text_hd.FontSize = 10;
end
P = polyfit(roll_return * 12, exce_return * 12, 1);
yfit = P(1) * (roll_return * 12) + P(2);
plot(roll_return * 12, yfit,'r-', 'linewidth', 1.5);
grid on;
xlabel('annulized roll return')
ylabel('annulized excess return')
set(gca, 'fontsize', 15)

% % q4
full_roll_returns = (-raw_data(:, 2:2:24) + raw_data(:, 1:2:23))./raw_data(:, 2:2:24);
cumed_return = zeros(size(returns, 1), 1);
temp_returns = full_roll_returns';
temp_full_returns = returns';

for i = 1:length(cumed_return)-1
    [~, I] = sortrows(temp_returns, i);
    end_four_idx = I(1:4);
    top_four_idx = I(end-3:end);
    cumed_return(i) = 0.25*sum((temp_full_returns(top_four_idx, i) - temp_full_returns(end_four_idx, i)));
end
figure
plot(cumprod(cumed_return+1), 'r-')
% close all
12 * mean(cumed_return)
std(cumed_return) * sqrt(12)
12 * mean(cumed_return)/(std(cumed_return) * sqrt(12))



% q5
temp_returns = returns + 1;
momen_return = zeros(size(temp_returns, 1) - 12, 1);
prod_mat = zeros(size(temp_returns, 1) - 12, 12);

for i = 13:size(temp_returns, 1)
    prod_mat(i-12, :) = prod(temp_returns(i-12:i-1, :));
    buffer_vector = prod_mat(i-12, :);

    [~, sort_max_index] = sort(buffer_vector, 'descend');
    max_index = sort_max_index(1:4);
    
    [~, sort_min_index] = sort(buffer_vector, 'ascend');
    min_index = sort_min_index(1:4);
    
    momen_return(i-12) = 0.25 * sum((temp_returns(i, max_index(:)) - temp_returns(i, min_index(:))));
end
cum_mat = cumprod(momen_return + 1);
figure
plot(cum_mat, 'r-')

12 * mean(momen_return)
std(momen_return) * sqrt(12)
12 * mean(momen_return)/(std(momen_return) * sqrt(12))