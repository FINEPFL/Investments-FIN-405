clearvars; close all; clc

% def helpers
getMean_A = @(x) mean(x);
getSD = @(x) std(x);
getSke = @(x) skewness(x);
getKurtosis = @(x) kurtosis(x);
getVaR95 = @(x) prctile(x, 5);

raw_data = xlsread('MSCI.xlsx');

% ---- structure of data table -----
% US | JP | CH | DE | Asia | Latin | EUME
px_table = raw_data(:, 2:3:20);

name_cell = cell(7, 1);
name_cell{1} = 'us';
name_cell{2} = 'jp';
name_cell{3} = 'ch';
name_cell{4} = 'de';
name_cell{5} = 'asia';
name_cell{6} = 'latin';
name_cell{7} = 'eu-me';

% compute monthly return of each region. note data is already monthly
return_table = diff(px_table)./px_table(1:(end-1), :);

%% (a)
for i=1:7
    disp('====================')
    disp(name_cell{i})
    sprintf('mean_A: %f', getMean_A(return_table(:, i)))
    sprintf('mean_G: %f', getMean_G(px_table(:, i)))
    sprintf('SD: %f',     getSD(return_table(:, i)))
    sprintf('Ske: %f',    getSke(return_table(:, i)))
    sprintf('Kur: %f',    getKurtosis(return_table(:, i)))
    sprintf('Var95: %f',  getVaR95(return_table(:, i)))
    sprintf('Min: %f',    min(return_table(:, i)))
    sprintf('Max: %f',    max(return_table(:, i)))
end

% hist(return_table(:, 1), 120)
% axis([-0.2, 0.4, 0 15])
figure()
subplot(2,4,1);
hist(return_table(:,1),100);
title('return US');
xlabel('return')
ylabel('frequency')
grid on


subplot(2,4,2)
hist(return_table(:,2),100);
title('return Japan');
xlabel('return')
ylabel('frequency')
grid on

subplot(2,4,3)
hist(return_table(:,3),100);
title('return Switzerland');
xlabel('return')
ylabel('frequency')
grid on

subplot(2,4,4)
hist(return_table(:,4),100);
title('return Germany');
xlabel('return')
ylabel('frequency')
grid on

subplot(2,4,5)
hist(return_table(:,5),100);
title('return Asia');
xlabel('return')
ylabel('frequency')
grid on

subplot(2,4,6)
hist(return_table(:,6),100);
title('return Latin');
xlabel('return')
ylabel('frequency')
grid on

subplot(2,4,7)
hist(return_table(:,7),100);
title('Returns EU--ME');
xlabel('return')
ylabel('frequency')
grid on

%% (b) verify rule of thumb M_G = M_A - 0.5*VAR
for i=1:7
    disp('====================')
    disp(name_cell{i})
    M_A = getMean_A(return_table(:, i));
    M_G = getMean_G(px_table(:, i));
    VAR = getSD(return_table(:, i))^2;
    sprintf('verify: %f', M_G-M_A+0.5*VAR)
end

%% (c)
clc
developted_matrix = return_table(:, 1:4);
COR_developted = corr(developted_matrix)

ME_matrix = return_table(:, 5:end);
COR_ME = corr(ME_matrix);

CH_ME_matrix = return_table(:, [3, 5, 6, 7]);
COR_CH_ME = corr(CH_ME_matrix)

%% (d)
win_size = 24;
list = [];

for i = 1:(size(return_table, 1)-24)+1
    cor_all(i) = get_avg_cor(return_table(i:(i-1)+win_size, :));
end
figure
plot(cor_all)
grid on;
xlabel('month index')
ylabel('averge correlation')
title('avg correlation in all markets, window\_size = 24')
set(gca, 'FontSize', 15)


