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

exampleData = csvread('example.csv')

%------------ PLOTTING FIGURE 2

N = 1e4; % number users to simulate

%required power at d0= 1m is capped at 20dB SNR_mac_maxdB




d= sqrt(rand(N,1)*R^2); %randomly created d - macro user

%{
t = 2*pi*rand(N,1);
r1 = d .* cos(t);
r2 = j*d .* sin(t);
d = sqrt(   abs(r1.^2) + abs(r2.^2)   );
%}


L = db2pow((sigma_sfdB)*randn(N,1)); %lognormal shadowing
%L = lognrnd(0,sigma_sfdB,[N,1]);

%P_txdB = SNR_maxdB + 10*pathloss_m*log10(d) + 10*log10(sigma_sfdB^2)


%P_txdB = -qfuncinv(alpha)   + ...
    %10*pathloss_m*log10(d);


    
shadowing_mean = (mean(d))^pathloss_m;   
%shadowing_mean = (R/2)^pathloss_m;
    
       
%P_tx = shadowing_mean * db2pow(( SNR_targetdB  - sigma_sfdB*qfuncinv(alpha) ));
P_tx = db2pow( SNR_targetdB -  sigma_sfdB*qfuncinv(alpha)    + 10*pathloss_m*log10(d)    )

%P_tx = db2pow(120);

mean(pow2db(P_tx))


%P_tx = db2pow(105);



P_rx = P_tx .* d.^(-pathloss_m) .* L; %macro, no wall loss, unity noise


%{
%impose mean rx received limit
change_meandB = mean(pow2db(P_rx)) - 20;
if change_meandB > 0
    P_rx = db2pow(pow2db(P_rx) - change_meandB);
end

%}

figure(1)

clf
hold on
cdfplot(pow2db(P_rx))
xlim([0 60]);
plot([0 60],[0.05 0.05], '--r')
plot([10 10],[0 1], '--r')
hold off


figure(2)
clf
h = histogram(pow2db(P_rx));
h.BinLimits = [0 60];





Pr_macro=P_rx;
bot_term = 1; %the noise power of the receiver, can be anything
SINR_macrolin= Pr_macro / bot_term;


mean(pow2db(SINR_macrolin))

capacity_macro = log(1 + SINR_macrolin);


x = [2 4 6 8]
y = [0.07 0.21 0.56 0.95]

figure(3)
clf
hold on
c = cdfplot(capacity_macro);
c.Color = 'k';
c.LineWidth = 2;
plot(x,y,'.r', 'MarkerSize', 10)


plot(exampleData(:,1), exampleData(:,2))

hold off

xlabel('C (bps/Hz)');
ylabel('cdf of C');
xlim([0 10]);
title('')



