%% David Dobbie
% 300340161
% Assignment 1, Question 8, ECEN 310


clear;
clf;
clc;

bandwidth=1;
cap = 0.0001:0.0001:4;

cdf_channel = 1 - exp(-(2.^(cap./bandwidth) - 1));



figure(1)

plot(cap, cdf_channel)
title({'Cumulative Distribution Function of Channel Capacity','with Exponentially Distributed SNR $$\mu=$$1'},'Interpreter','latex')
xlabel('Channel capacity per bandwidth $$\frac{C}{B}$$ (bits)','Interpreter','latex');
ylabel('Capacity Guaranteed Probability','Interpreter','latex');




