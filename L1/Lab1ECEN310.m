% David Dobbie : 300340161
% ECEN 310 / ENGR 440 Communications Engineering
% Lab 1 - Bit Error Rate

clear all; close all; clc;

%Q1:

set(0, 'defaulttextInterpreter','latex')

% init MPSK
M = 2;
constel = exp((j * 2* pi * (0:M-1))/ (M))

N = 1e3;
num_tests= 50;


s = zeros(N,1);
n = zeros(N,1);
r = zeros(N,1);
sest = zeros(N,1);

Es = 1;

results = zeros(num_tests,2);
iter = 1;

figure(1)

SNRdB_axis = 0:5:15;

for SNRdBindx = 1:length(SNRdB_axis);
    SNRdB = SNRdB_axis(SNRdBindx);
    No = Es/db2pow(SNRdB);

    rng(6) % sets generator seed

    s = randsrc(N,1,constel); % get rnd symbols, tx
    n = sqrt(No/2)*complex(randn(N,1),randn(N,1)); % noise samples
    r = s + n; % rx


    for indx = 1:N
        %returns decision point closest to the received message
        [dmin, const_indx] = min(abs(r(indx) - constel)); 
        sest(indx) = constel(const_indx);
    end

    subplot(2,2,SNRdBindx)
    
    hold on
    plot(real(constel)', imag(constel)', 'r+', 'linewidth', 2, 'markersize' , 8)
    plot(real(r)', imag(r)', 'b.')
    hold off
    
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis([-2 2 -2 2])
    grid on
    
    xlabel('$ \psi_1 (t)$')
    ylabel('$ \psi_2 (t)$')
    
    title(['$SNR_{dB}$ = ' num2str(SNRdB) 'dB'])

    
    %[s r sest];

    %SERnumlog = log10(nnz(s-sest)) - log10(N); %computation in log to deal with small num
    %SERnum = 10^SERnumlog; %convert to linear scale
    SERnum = nnz(s-sest)/N;
    BERnum = SERnum *log2(M); %since 4-ary system (2 bits per symbol)
    
    
end



