clearvars; close all; clc
raw_data = load('raw_data.csv');

coupon_rate   = raw_data(:, 1);
maturity_year = raw_data(:, 2);
current_price = raw_data(:, 3);
par_value     = raw_data(:, 4);

% getting the payment matrix
N = size(coupon_rate, 1);
for i = 1:N
    annual = 2 * maturity_year(i); % data in half year
    Y(i, 1:annual) = coupon_rate(i)/100 * par_value(i)/2;
    Y(i, annual) = par_value(i) + Y(i, annual);
end

% according to formula on slide 33, we can get discount factor
DF = Y \ current_price;
bond_price = DF .* par_value;
for i = 1:N
    ZCR(i) = 2 * ((100/bond_price(i))^(0.5/maturity_year(i)) - 1);
end
figure
plot(maturity_year, ZCR, 'b','linewidth', 3); grid on;
title('ZCR curve')
set(gca, 'fontsize', 12)
disp(ZCR')
    


