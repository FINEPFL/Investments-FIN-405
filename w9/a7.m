clearvars; close all; clc

raw_data = xlsread('BAML_total_return.xlsx');

rf    = raw_data(:, end);
block = raw_data(:, 2:end-1);
N     = size(block, 1);
T     = size(block, 2);

% annulized = zeros(N-1, T);
% for i = 2:N
%     annulized(i-1, :) = ((block(i, :) ./ block(i-1, :)) .^ 12 - 1) .* 100;
% end
% to compare
% A = ((annulized(:,1:4)./100+1).^(1/12) - 1)*100;

returns = (block(2:end, :) ./ block(1:end-1, :) - 1) .* 100;
returns = [returns rf(1:end-1)];

avg_exc_ret = mean(returns*12)./100 - mean(returns(:, 5))/100
Sigma       = sqrt(12) * std(returns) ./ 100
SR          = avg_exc_ret ./ Sigma