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




capacity_macro = generateCM(100, 0, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto5 =  generateCM(100, 5, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto100 =  generateCM(100, 100, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);
capacityM_femto300 =  generateCM(100, 300, R, R_o, r_f, sigma_sfdB, SNR_targetdB, SNR_maxdB,  pathloss_m, NodB, WdB);


capacityF_femto5 =  generateFM(10, 5, R, R_o, r_f, sigma_sfdB, SNR_targetdB, ...
         SNR_maxdB,  pathloss_m, NodB, WdB, p_act, p_ind);


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

c = cdfplot(capacityF_femto5);
c.Color = 'b';
c.LineStyle = '-';
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
%% MACRO USER CAPACITY
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
    
    
    % generate femtos range
    magnitude_femto = sqrt(abs(rand(N_femto,1)*R_outRange^2));
    bearing_femto = 2*pi*(rand(N_femto,1));
    pos_femto = magnitude_femto .* exp(1i*bearing_femto);

    
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
    P_tx = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.975)  + 10*pathloss_m*log10(mean(dm)) + NodB) %for macro signal strength
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
                + 10*pathloss_m*log10((2/3)*r_femto) + NodB); %compute mean dist as 2/3 radius
            L_femto = db2pow((sigma_sfdB)*randn(N,1)); %lognormal shadowing
            P_rx_femto = (P_tx_femto .* dist_from_femto.^(-pathloss_m) .* L_femto) / db2pow(wallLossdB);
            
            % maxes out received power from the femto for each user to 20dB
            for idx = 1:length(P_rx_femto)
                if pow2db(P_rx_femto(idx)) - NodB > SNR_maxdB
                    P_rx_femto(idx) = db2pow(SNR_maxdB + NodB) ;
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








%% FEMTO USER CAPACITY
function [capacity] = generateFM(macro_userDens, femto_Dens, R , R_outRange, r_femto, ...
    sigma_sfdB, SNR_targetdB, SNR_maxdB, pathloss_m ,NodB, wallLossdB, p_activity , p_indoor)
    lambda_macro = macro_userDens * pi * (R*1e-3)^2; % convert to km
    N = poissrnd(lambda_macro); % number users to simulate

    
    % generate macro average to set macro transmit power
    magnitude = sqrt(abs(rand(N,1)*R^2));
    bearing = 2*pi*(rand(N,1));
    pos = magnitude .* exp(1i*bearing);


    lambda_femto = femto_Dens * pi * (R_outRange*1e-3)^2; % convert to km
    N_femto = poissrnd(lambda_femto); % number users to simulate
    
    % generate femtos
    magnitude_femto = sqrt(abs(rand(N_femto,1)*R_outRange^2));
    bearing_femto = 2*pi*(rand(N_femto,1));
    pos_femto = magnitude_femto .* exp(1i*bearing_femto);
   
    %{
    % generate femto users
    magnitude_femto_user = sqrt(abs(rand(N_femto,1)*r_femto^2));
    bearing_femto_user = 2*pi*(rand(N_femto,1));
    pos_fromFC_femto_user = magnitude_femto_user .* exp(1i*bearing_femto_user);
    
    pos_femto_user = pos_fromFC_femto_user + pos_femto; %adds in position relative to macro (our origin)
    
    magnitude_FU_from_macro = abs(pos_femto_user);
    
    isInside = [1 0];
    isFUInside = randsrc(N_femto ,1, isInside); % is the FU inside 0.5 prob    
    %}
    dm = magnitude;  %from macro to macro user

    


    SINR_femto_user = zeros(N_femto,1);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    I = 0; %femto macro interference
    % add in femto interference to the macro user
%{
      for FU_idx = 1:N_femto % for each femto user
            FU_own_cell_P_rx = 0;

            magnitude_femto_user = sqrt(abs(rand(N_femto,1)*r_femto^2));
            bearing_femto_user = 2*pi*(rand(N_femto,1));
            pos_fromFC_femto_user = magnitude_femto_user .* exp(1i*bearing_femto_user);

            pos_femto_user = pos_fromFC_femto_user + pos_femto(; %adds in position relative to macro (our origin)

            magnitude_FU_from_macro = abs(pos_femto_user);

            isInside = [1 0];
            isFUInside = randsrc(N_femto ,1, isInside); % is the FU inside 0.5 prob   


            L = db2pow((sigma_sfdB)*randn(N_femto,1)); %lognormal shadowing

            %------------ SETS TRANSMIT POWER ACCORDING TO eq 2
            %sets transmit power macro such that mean macro SNR  meets threshold (eq 2)
            P_tx_m = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.95)  + 10*pathloss_m*log10(mean(dm)) + NodB) 
            % caps received power for femto users
            P_rx_m_to_fu = (P_tx_m .* magnitude_FU_from_macro.^(-pathloss_m) .* L);

            % cap the macro transmit power of the macro to the femto user
            for idx = 1:length(P_rx_m_to_fu)
                if pow2db(P_rx_m_to_fu(idx)) - NodB > SNR_maxdB
                    P_rx_m_to_fu(idx) = db2pow(SNR_maxdB + NodB) ;
                end 
            end


                
              
            for idx_femto = 1:N_femto %for each femto
                
 
    %{
                figure(8)
                clf
                hold on
                plot(real(pos_femto), imag(pos_femto), '+r')
                plot(real(pos_femto_user), imag(pos_femto_user), '.g')
                plot(0,0, '+k', 'MarkerSize', 20)
                hold off
      %}          
                
                pos_fromFC_femto_user = abs( pos_femto(idx_femto) - pos_femto_user(FU_idx));   % distance between femto user and the femto cell

                % setting femto transmit power to attempt 95% connectivity within its
                % range for received users
                P_tx_femto_cell = db2pow(+SNR_targetdB    -   sigma_sfdB*qfuncinv(0.95) ...
                    + 10*pathloss_m*log10((2/3)*r_femto) + NodB); %compute mean dist as 2/3 radius

                L_femto = db2pow((sigma_sfdB)*randn(1,1)); %lognormal shadowing

                P_rx_femto = (P_tx_femto_cell .* pos_fromFC_femto_user.^(-pathloss_m) ...
                    .* L_femto) ./ db2pow(wallLossdB.*(isFUInside(FU_idx)+1));

                %pow2db(P_rx_femto);
                
                % maxes out received power for each user to 20dB
                for idx = 1:length(P_rx_femto)
                    if pow2db(P_rx_femto(idx)) - NodB > SNR_maxdB
                        P_rx_femto(idx) = db2pow(SNR_maxdB + NodB) ;
                    end 
                end
                
                figure(99)
                clf
                cdfplot(pow2db(P_rx_m_to_fu) - NodB)
                

                rayleigh = (  abs(    sqrt(1/2) * ( randn(1) + 1j*randn(1) )    )  ).^2; %unit variance, zero mean
                % creates interference term with the received power from the
                % femto
                interfere_femto = (P_rx_femto .* rayleigh);

                if(FU_idx == idx_femto) % if it is the femto that the FU is using, it is NOT interference
                    FU_own_cell_P_rx = interfere_femto; % the femto users own femto
                    interfere_femto = 0; 
                end

                I = I + interfere_femto;

            end
            rayleigh = (  abs(    sqrt(1/2) * ( randn(1) + 1j*randn(1) )    )  ).^2; %unit variance, zero mean
            
            macro_interfere = P_rx_m_to_fu(FU_idx) * abs(pos_femto_user(FU_idx)) * rayleigh *...
                L(FU_idx) / db2pow(wallLossdB.*isFUInside(FU_idx));
            macro_interfere = 0;
            
            SINR_femto_user(FU_idx) = FU_own_cell_P_rx./(db2pow(NodB) + I + macro_interfere);
      end 
    %}
  
  
  
    capacity = log2(1 + SINR_femto_user);
end

















