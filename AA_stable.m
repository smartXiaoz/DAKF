close all;clear all; clc
N = 1000;     %total number of time steps
n = 2;
m = 4;
F = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)]; % 模型
H = ones(m,n);
epsilon = 1e-6;
del_T=0.1;

xe_ini =ones(n,1)*4; %估计初始值
xe_MEE = xe_ini;




Iter_num = 10;