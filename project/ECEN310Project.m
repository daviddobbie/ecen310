%% ECEN 310 Project
% David Dobbie - 300340161

clc
clear all
set(0,'defaultTextInterpreter','latex');

r_f = 20;
R = 1000;
R_o = 2e3;
p_act = 0.25;
p_ind = 0.5;
pathloss_m = 3;
pathloss_f = 3;
sigma_sfdB =8;
SNR_targetdB = 10;
alpha = 0.95;
SNR_maxdB = 20;
WdB = 10;
NodB = 130;

exampleData = csvread('example.csv');

%------------ PLOTTING FIGURE 2

N = 1e5; % number users to simulate





%required power at d0= 1m is capped at 20dB SNR_mac_maxdB

d= sqrt(abs(rand(N,1)*R^2)); %randomly created d - macro user
L = db2pow((sigma_sfdB)*normrnd(0,1,N,1)); %lognormal shadowing


%------------ SETS TRANSMIT POWER ACCORDING TO eq 2
%sets transmit power such that mean macro SNR  meets threshold (eq 2)
P_tx = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.975)  + 10*pathloss_m*log10(mean(d)) + NodB) 


P_rx = (P_tx .* d.^(-pathloss_m) .* L);


for idx = 1:length(P_rx)
    if pow2db(P_rx(idx)) - NodB > SNR_maxdB
        P_rx(idx) = db2pow(SNR_maxdB + NodB) ;
    end 
end




figure(11)
clf
SNR_rx_dB = pow2db(P_rx) - NodB;
cdfplot(SNR_rx_dB)
xlabel('Mean $SNR_m$')
ylabel('cdf of mean $SNR_m$')
title('')

I = 0; %femto macro interference

%rayleigh distribution on receiver
rayleigh = (  abs(    sqrt(1/2) * ( normrnd(0,1,N,1) + 1j*normrnd(0,1,N,1) )    )  ).^2; %unit variance, zero mean
%rayleigh =1;
SINR_macro_dB = pow2db(P_rx) + pow2db(rayleigh) - pow2db((db2pow(NodB) + I));
SINR_macro = db2pow(SINR_macro_dB);
%SINR_macro = P_rx.* rayleigh / (db2pow(NodB) + I); %instant rx power macro only


figure(1)

clf
hold on
cdfplot(pow2db(SINR_macro))
xlim([0 60]);
plot([0 60],[0.05 0.05], '--r')
plot([10 10],[0 1], '--r')
hold off
xlabel('Instantenous for each macro user $SINR_m$')
ylabel('cdf of $SINR_m$')

capacity_macro = log2(1 + SINR_macro);


figure(3)
clf
hold on
c = cdfplot(capacity_macro);
c.Color = 'k';
c.LineWidth = 2;
plot(exampleData(:,1), exampleData(:,2),'LineWidth', 1.5)

hold off

xlabel('C (bps/Hz)');
ylabel('cdf of C');
xlim([0 10]);
title('')


figure(1)
hold on
plot(pow2db(2.^exampleData(:,1) - 1), exampleData(:,2))
hold off


