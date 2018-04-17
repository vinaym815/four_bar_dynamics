% close all;
clear variables;
% clc;
%%
global thStart t0 tf
L0 = 3;  L1 = 9;  L2 = 10;  L3 = 7;
m1 = 1; m2 = 1; m3 = 1;
tau = 6;

th0 = 1.5708;
t0 = 0;     tf = 15;

thStart = th0;

sim('fourbarConstTau.slx')
plot(theta,'LineWidth',2)