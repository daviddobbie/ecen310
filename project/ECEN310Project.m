%% ECEN 310 Project
% David Dobbie - 300340161

clc
clear all
set(0,'defaultTextInterpreter','latex');

r_f = 20;
R = 1e3;
R_o = 2e3;
p_act = 0.25;
p_ind = 0.5;
pathloss_m = 3;
pathloss_f = 3;
sigma_sfdB = 8;
SNR_targetdB = 10;
alpha = 0.95;
SNR_maxdB = 20;
WdB = 10;



%------------ PLOTTING FIGURE 2

N = 1e6; % number users to simulate

%required power at d0= 1m is capped at 20dB SNR_mac_maxdB




d= sqrt(rand(N,1)*R^2); %randomly created d - macro user
XdB = sigma_sfdB*randn(N,1); %lognormal shadowing


%P_txdB = SNR_maxdB + 10*pathloss_m*log10(d) + 10*log10(sigma_sfdB^2)


%P_txdB = -qfuncinv(alpha)   + ...
    10*pathloss_m*log10(d);

P_txdB = 109


SNR_rxdB = P_txdB + XdB - 10*pathloss_m*log10(d);

if SNR_rxdB > SNR_maxdB
    SNR_rxdB = SNR_maxdB;
end


figure(1)

clf
hold on
cdfplot(SNR_rxdB)
xlim([0 60]);
plot([0 60],[0.05 0.05], '--r')
plot([10 10],[0 1], '--r')
hold off








Pr_macro=SNR_rxdB;
bot_term = sigma_sfdB^2;
SINR_macrolin= 10.^(Pr_macro./10) / bot_term;


capacity_macro = log2(1 + SINR_macrolin);



figure(2)
clf
hold on
cdfplot(capacity_macro)

hold off

xlabel('C (bps/Hz)');
ylabel('cdf of C');
xlim([0 10]);

