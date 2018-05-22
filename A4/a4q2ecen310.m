
clear all
close all
clf
clc 

SNRdB = 0:1:10;

SNRlin = 10.^(SNRdB./10)


BER_orig = 2*qfunc(sqrt(2* SNRlin))

BER_new = 1 - sqrt(SNRlin./(1+SNRlin))


figure(1)
hold on
loglog(SNRdB, BER_orig)
loglog(SNRdB, BER_new)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off


