%*************** Q3 and Q4 for investments assignment 2 ******************
% Authors:
%     Mengjie Zhao,
%     Tianxiao Ma,
%     Wenfei Chu
%
% Contacts:
%     {mengjie.zhao, tianxiao.ma, wenfei.chu}@epfl.ch
%
% Please feel free to drop emails to the authors if you met any problem in
% evaluating the codes.

clearvars; close all; clc

%% Q3 Alternative portfolio strategies

% processing the raw data
% Treasurey | Corp Bonds | US stock | Commodies | 1-month T bill
raw_data  = xlsread('totalReturns.xlsx');
price_mat = raw_data(:, 1:4);

% calculate return from price data
R  = diff(price_mat)./price_mat(1:end-1, :);
rf = raw_data(1:end, 5); % the risk free rate

% def parameters of sliding  
% 10 years window size
win_size = 120;

% As the data is monthly return and the portfolio will be hold for one
% month, then the step will be one. Hence,  (size(R, 1) - win_size) steps
steps = length(R) - win_size;

% Initialize the memory to make it run faster
weight_TAN = zeros(steps, 4);
weight_GMV = zeros(steps, 4);
weight_RP  = zeros(steps, 4);
weight_EW  = 0.25 .* ones(steps, 4); 
port_R     = zeros(steps, 4);

rolling_matrix = zeros(win_size, 4);

% def getters
getSigma = @(return_matrix) cov(return_matrix);
getMu    = @(return_matrix) mean(return_matrix)'; % avg over vertical axis
getA     = @(Sigma) ones(1, 4) * inv(Sigma) * ones(4, 1);
getB     = @(Sigma, Mu) ones(1, 4) * inv(Sigma) * Mu; % Note the input Mu 
                                                      % is a already 
                                                      % column vector.
getC     = @(Sigma, Mu) Mu' * inv(Sigma) * Mu;
getDelta = @(A, B, C) A*C - B^2;
                       
for i = 1:steps
    
    rolling_matrix = R(i:i+win_size-1, :);
    
    mu    = getMu(rolling_matrix);
    sigma = getSigma(rolling_matrix);
    A     = getA(sigma);
    B     = getB(sigma, mu);
    C     = getC(sigma, mu);
    Delta = getDelta(A, B, C);
    rf_   = rf(i+win_size+1);
    
    weight_TAN(i, :) = (inv(sigma) * (mu - rf_/12 * ones(4, 1))) /...
                                                    (B - A * rf_/12);

    weight_GMV(i, :) = inv(sigma) * ones(4, 1) ./ ...
                                    (ones(1, 4) * inv(sigma) * ones(4, 1));
                                
    stddev = std(rolling_matrix);
    weight_RP(i, :) = (1./stddev)/sum(1./stddev);
    
    port_R(i, 1) = weight_TAN(i, :) * R(i+(win_size), :)';
    port_R(i, 2) = weight_GMV(i, :) * R(i+(win_size), :)';
    
    port_R(i, 3) = weight_RP(i, :) * R(i+(win_size), :)';
    port_R(i, 4) = weight_EW(1, :) * R(i+(win_size), :)';
 
end

figure(1)
plot(1:steps, weight_TAN, 'linewidth', 1.5)
xlabel('monthly index', 'interpreter', 'latex')
ylabel('weight', 'interpreter', 'latex')
title('TAN', 'interpreter', 'latex')
legend('Treasurey', 'Corp Bonds', 'US stock', 'Commodies');
set(gca, 'fontsize', 15)

figure(2)
plot(1:steps, weight_GMV, 'linewidth', 1.5)
xlabel('monthly index', 'interpreter', 'latex')
ylabel('weight', 'interpreter', 'latex')
title('GMV', 'interpreter', 'latex')
legend('Treasurey', 'Corp Bonds', 'US stock', 'Commodies');
set(gca, 'fontsize', 15)

figure(3)
plot(1:steps, weight_RP, 'linewidth', 1.5)
xlabel('monthly index', 'interpreter', 'latex')
ylabel('weight', 'interpreter', 'latex')
title('RP', 'interpreter', 'latex')
legend('Treasurey', 'Corp Bonds', 'US stock', 'Commodies');
set(gca, 'fontsize', 15)

figure(4)
plot(1:steps, weight_EW, 'linewidth', 1.5)
xlabel('monthly index', 'interpreter', 'latex')
ylabel('weight', 'interpreter', 'latex')
legend('Treasurey', 'Corp Bonds', 'US stock', 'Commodies');
title('EW', 'interpreter', 'latex')
set(gca, 'fontsize', 15)

% b) compute mean and std of portfolio returns port_R, and sharp ratio
port_mean = mean(port_R)
port_std = std(port_R)
rf_mean = mean(rf(win_size+1:end) / 12);
sharp_ratio = (port_mean - rf_mean) ./ port_std .* sqrt(12)

% c) explanation

% d) Plot the minimum-variance frontier, given the mean of portfolio, we can
% calculate the variance using the mean-variance relation (Au^2-2Bmu+C)/Dt
d_price_mat  = raw_data(win_size+1:end, 1:4);
d_return_mat = diff(d_price_mat)./d_price_mat(1:end-1, :);

d_mu    = getMu(d_return_mat);
d_sigma = getSigma(d_return_mat);
d_A     = getA(d_sigma);
d_B     = getB(d_sigma, d_mu);
d_C     = getC(d_sigma, d_mu);
d_Delta = getDelta(d_A, d_B, d_C);

mu_axis  = 0.001:0.00001:0.012;
var_axis = sqrt((d_A * mu_axis.^2 - 2 * d_B * mu_axis + d_C)./d_Delta);
figure(5)
plot(var_axis, mu_axis); grid on; hold on;
set(gca, 'fontsize', 15)
xlabel('Standard Deviation $\sigma$', 'interpreter', 'latex')
ylabel('Expected Return $\mu$', 'interpreter', 'latex')

tan_mean = 0.0052; gmv_mean = 0.0055; rp_mean = 0.0061; ew_mean = 0.0060;
tan_std  = 0.0176; gmv_std  = 0.0125; rp_std  = 0.0137; ew_std  = 0.0213;

plot(tan_std, tan_mean,'*', 'markersize', 5)
plot(gmv_std, gmv_mean,'^', 'markersize', 5)
plot(rp_std, rp_mean,'o', 'markersize', 5)
plot(ew_std, ew_mean,'s', 'markersize', 5)
legend('mean-std locus', 'TAN', 'GMV', 'RP', 'EW')

% e) 1 dollar performance of all strategies
perform = -1 .* ones(steps, 5);
for i=1:steps
    if i == 1
        perform(i, 1:4) = 1 + port_R(i, :);
        perform(i, 5)   = 1 + rf(i+120)/12;
    else
        perform(i, 1:4) = perform(i-1, 1:4) .* (1 + port_R(i, :));
        perform(i, 5)   = perform(i-1, 5) .* (1 + rf(i+120)/12);
    end
end
figure(6)
plot(1:steps, perform, 'linewidth', 1.5)
xlabel('monthly index', 'interpreter', 'latex')
ylabel('Value(\$)')
title('1\$ cumulative performance', 'interpreter', 'latex')
legend('TAN', 'GMV', 'RP', 'EW', 'rf')
set(gca, 'fontsize', 15)        
   