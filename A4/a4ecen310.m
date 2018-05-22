%% ECEN 310 Assignment 4
% David Dobbie 300340161


%Q1)
% plotting SNR gain
clear all;
close all;
clf;
clc


L = 1:8;



set(0,'defaulttextInterpreter','latex') 


gainEGC = (1 + (L-1)*(pi/4))' *100

gainSC = ones(length(L),1);
for el = L
    val = 1:el;
    gainSC(el) = sum(1./val)' * 100;
end
gainSC

gainMRC = L' *100


figure(1)
hold on
plot(L,gainEGC,'LineWidth',2,'LineStyle', '--');
plot(L,gainSC,'LineWidth',2,'LineStyle', '--');
plot(L,gainMRC,'LineWidth',2,'LineStyle', '--');
title('Mean SNR gain for each combining scheme')
xlabel('Antenna Count L')
ylabel('Percentage of gain compared to $\bar{\gamma}$ [\%]')

legend('Equal Gain Combining','Selection Combining','Mixed Ratio Combining')




