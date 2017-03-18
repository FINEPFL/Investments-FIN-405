clearvars; close all; clc

%% Q1
mu    = [0.1, 0.12, 0.11]';
R_std = [0.1 ,0.11, 0.15];
R_var = [0.1 ,0.11, 0.15]' * [0.1 ,0.11, 0.15];

B  = [0.2, 0.5,   0;
      0.4, 0.6,   0;
      0  , 0  , 0.3;];
  
F_std = [0.1, 0.15, 0.2]';
F_var = diag(F_std.^2);

product  = B' * F_var * B;
epsl_var = diag(diag(R_var)-diag(product))
% a)
sqrt(product)
sqrt(epsl_var)

% b)
rho1_2 = product(1, 2)/(R_std(1) * R_std(2))

% c)
w = [0.25, 0.45, 0.3]';
w' * B'

% d) Denote I as the diag matrix, meaning exposure for that factor is 1
% Then we have B * W = I
I = eye(3);
W = inv(B) * I

%% Q2
clc
% a)
inv([1, 1, 0.5;
     1, 3, 0.2;
     1, 3, -0.5;]) * [0.12;
                      0.134;
                      0.12;]
                  
% b) skiped. Pls refer to docs      
