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

%{


%------------ PLOTTING FIGURE 2

mean_macro_user_density = 1000; %users per km^2
lambda_macro = mean_macro_user_density * pi * (R*1e-3)^2 % convert to km
N = poissrnd(lambda_macro) % number users to simulate

%required power at d0= 1m is capped at 20dB SNR_mac_maxdB
% generate macro users
% randomly intialised polar coordinates
magnitude = sqrt(abs(rand(N,1)*R^2));
bearing = 2*pi*(rand(N,1));
pos = magnitude .* exp(1i*bearing);

d= magnitude; %randomly created d - macro user
L = db2pow((sigma_sfdB)*randn(N,1)); %lognormal shadowing

%------------ SETS TRANSMIT POWER ACCORDING TO eq 2
%sets transmit power such that mean macro SNR  meets threshold (eq 2)
P_tx = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.975)  + 10*pathloss_m*log10(mean(d)) + NodB) 
P_rx = (P_tx .* d.^(-pathloss_m) .* L);
for idx = 1:length(P_rx)
    if pow2db(P_rx(idx)) - NodB > SNR_maxdB
        P_rx(idx) = db2pow(SNR_maxdB + NodB) ;
    end 
end

I = 0; %femto macro interference
%rayleigh distribution on receiver
rayleigh = (  abs(    sqrt(1/2) * ( randn(N,1) + 1j*randn(N,1) )    )  ).^2; %unit variance, zero mean
SINR_macro_dB = pow2db(P_rx) + pow2db(rayleigh) - pow2db((db2pow(NodB) + I));
SINR_macro = db2pow(SINR_macro_dB);

capacity_macro = log2(1 + SINR_macro);

%}

capacity_macro = generateCM(1000, 0, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto5 =  generateCM(1000, 5, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto100 =  generateCM(1000, 100, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto300 =  generateCM(1000, 300, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);

figure(2)
clf
hold on
c = cdfplot(capacity_macro);
c.Color = 'k';
c.LineWidth = 2;

c = cdfplot(capacityM_femto5);
c.Color = 'r';
c.LineWidth = 2;

c = cdfplot(capacityM_femto100);
c.Color = 'r';
c.LineStyle = '--';
c.LineWidth = 2;

c = cdfplot(capacityM_femto300);
c.Color = 'r';
c.LineStyle = ':';
c.LineWidth = 2;

hold off
xlabel('C (bps/Hz)');
ylabel('cdf of C');
xlim([0 10]);
title('')


%% functions
% inputs:
% results: we find the capacity of each macro user due to macro and femto
% inteference
function [capacity] = generateCM(macro_userDens, femto_Dens, R , R_outRange, r_femto, ...
    sigma_sfdB, SNR_targetdB, SNR_maxdB, pathloss_m ,NodB, wallLossdB)
    lambda_macro = macro_userDens * pi * (R*1e-3)^2; % convert to km
    N = poissrnd(lambda_macro); % number users to simulate
    
    lambda_femto = femto_Dens * pi * (R*1e-3)^2; % convert to km
    N_femto = poissrnd(lambda_femto); % number users to simulate

    %required power at d0= 1m is capped at 20dB SNR_mac_maxdB
    % generate macro users
    % randomly intialised polar coordinates
    magnitude = sqrt(abs(rand(N,1)*R^2));
    bearing = 2*pi*(rand(N,1));
    pos = magnitude .* exp(1i*bearing);

    lambda_femto = femto_Dens * pi * (R_outRange*1e-3)^2; % convert to km
    N_femto = poissrnd(lambda_femto); % number users to simulate
    
    
    % generate femtos inside range
    magnitude_femto = sqrt(abs(rand(N_femto,1)*R_outRange^2));
    bearing_femto = 2*pi*(rand(N_femto,1));
    pos_femto = magnitude_femto .* exp(1i*bearing_femto);
    
    % generate femtos outside range
    
    %{
    magnitude_femto_outside = sqrt(abs(rand(N_femto,1)*R^2)) + (R_outRange - R);
    bearing_femto_outside = 2*pi*(rand(N_femto,1));
    pos_femto_outside = magnitude_femto_outside .* exp(1i*bearing_femto_outside);
    
    pos_femto = [pos_femto ;  pos_femto_outside ]
%}
    
    figure(8)
    clf
    hold on
    plot(real(pos), imag(pos), '.b')
    plot(real(pos_femto), imag(pos_femto), '+r')
    hold off
    
    dm= magnitude; %randomly created d - macro user
    L = db2pow((sigma_sfdB)*randn(N,1)); %lognormal shadowing

    %------------ SETS TRANSMIT POWER ACCORDING TO eq 2
    %sets transmit power such that mean macro SNR  meets threshold (eq 2)
    P_tx = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.975)  + 10*pathloss_m*log10(mean(dm)) + NodB) 
    P_rx = (P_tx .* dm.^(-pathloss_m) .* L);
    for idx = 1:length(P_rx)
        if pow2db(P_rx(idx)) - NodB > SNR_maxdB
            P_rx(idx) = db2pow(SNR_maxdB + NodB) ;
        end 
    end

    I = 0; %femto macro interference
    
    % add in femto interference to the macro user
    if N_femto > 0
      for idx = 1:N_femto
          dist_from_femto = zeros(N,1);
          for user = 1:N %distance between each macro user and the current femto
            dist_from_femto(user) = abs( pos_femto(idx) - pos(user));                
          end
            % setting femto transmit power to attempt 95% connectivity within its
            % range for received users
            P_tx_femto = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.95) ...
                + 10*pathloss_m*log10((2/3)*r_femto) + NodB);
            L_femto = db2pow((sigma_sfdB)*randn(N,1)); %lognormal shadowing
            P_rx_femto = (P_tx_femto .* dist_from_femto.^(-pathloss_m) .* L_femto) / wallLossdB;
            
            % maxes out received power for each user to 20dB
            for idx = 1:length(P_rx)
                if pow2db(P_rx(idx)) - NodB > SNR_maxdB
                    P_rx(idx) = db2pow(SNR_maxdB + NodB) ;
                end 
            end
            
            rayleigh = (  abs(    sqrt(1/2) * ( randn(N,1) + 1j*randn(N,1) )    )  ).^2; %unit variance, zero mean
            % creates interference term with the received power from the
            % femto
            interfere_femto = (P_rx_femto .* rayleigh);
            I = I + interfere_femto;
        end 
    end
    
    
    
    
    %rayleigh distribution on receiver
    rayleigh = (  abs(    sqrt(1/2) * ( randn(N,1) + 1j*randn(N,1) )    )  ).^2; %unit variance, zero mean
    SINR_macro_dB = pow2db(P_rx) + pow2db(rayleigh) - pow2db((db2pow(NodB) + I));
    SINR_macro = db2pow(SINR_macro_dB);

    capacity = log2(1 + SINR_macro);
end