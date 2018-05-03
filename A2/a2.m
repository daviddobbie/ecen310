%% David Dobbie
% 300340161
% Assignment 2, Question 4, ECEN 310

close all;
clear all;
clc;

set(0, 'defaulttextInterpreter','latex')


constel = [-3 -1 1 3]  % init signal constell
M = 4;

N = 1e7;



num_tests= 50;


s = zeros(N,1);
n = zeros(N,1);
r = zeros(N,1);
sest = zeros(N,1);

Es = 1;

results = zeros(num_tests,2);
iter = 1;



for SNRdB = logspace(0,2,num_tests);
    No = Es/db2pow(SNRdB);

    rng(6) % sets generator seed

    s = randsrc(N,1,constel); % get rnd symbols, tx
    n = sqrt(No/2)*randn(N,1); % noise samples
    r = s + n; % rx


    for indx = 1:N
        %returns decision point closest to the received message
        [dmin, const_indx] = min(abs(r(indx) - constel)); 
        sest(indx) = constel(const_indx);
    end


    %[s r sest];

    %SERnumlog = log10(nnz(s-sest)) - log10(N); %computation in log to deal with small num
    %SERnum = 10^SERnumlog; %convert to linear scale
    SERnum = nnz(s-sest)/N;
    BERnum = SERnum *2; %since 4-ary system (2 bits per symbol)
    results(iter, :) = [SNRdB BERnum];
    
    iter = iter + 1;
end


figure(1)

loglog(results(:,1), results(:,2),'LineWidth',6)
xlabel('SNR/bit (dB/bit)');
ylabel('Bit Error Rate');
title('Comparsion of Experimental and Theoretical BER for a 4-ary PAM Constellation')
hold on
axis([1e0 10^(1.2) 1e-6 1e0])
grid on


% plots the analytical expression found in class
iter = 1;
analytData = zeros(num_tests,2);
for SNRdB = logspace(0,2,num_tests);
    No = Es/db2pow(SNRdB);
    analytSER = 2*((M-1)/M)* qfunc(sqrt((2*Es)/No));
    
    analytData(iter, :) = [SNRdB analytSER];
    
    iter = iter + 1;
end
%analytBER = 2*analytSER;
analytData(:,2) = 2*analytData(:,2); %convert to BER

loglog(analytData(:,1), analytData(:,2),'LineWidth',3);

% plots the union bound accounding to the min distance approx

mindist = Es *2; %due to M-ary PAM around centre
iter = 1;
unionBoundQFunc = zeros(num_tests,2);
for SNRdB = logspace(0,2,num_tests)
    No = Es/db2pow(SNRdB);
    unionQfuncSER = (M-1)* qfunc(mindist/sqrt(2*No)); 
    unionQfuncBER = 2*unionQfuncSER;
    unionBoundQFunc(iter, :) = [SNRdB unionQfuncSER];
    iter = iter + 1;
end
unionBoundQFunc(:,2) = 2*unionBoundQFunc(:,2); %convert to BER
loglog(unionBoundQFunc(:,1), unionBoundQFunc(:,2),'LineWidth',3);



% plots the nbound using min dist approx and bound on q func

mindist = Es *2; %due to M-ary PAM around centre
iter = 1;
unionBoundQFuncApprox = zeros(num_tests,2);
for SNRdB = logspace(0,2,num_tests)
    No = Es/db2pow(SNRdB);
    unionQfuncSER = ((M-1)/2) * exp(-(mindist^2)/(4*No)) ;
    unionQfuncBER = 2*unionQfuncSER;
    unionBoundQFuncApprox(iter, :) = [SNRdB unionQfuncSER];
    iter = iter + 1;
end
unionBoundQFuncApprox(:,2) = 2*unionBoundQFuncApprox(:,2); %convert to BER
loglog(unionBoundQFuncApprox(:,1), unionBoundQFuncApprox(:,2),'LineWidth',3);

hold off
lgnd = legend('Simulated $ P_m $' , 'Analytical Expression', 'Union bound min dist.', 'Union bound min dist. approx')
lgnd.Location = 'southwest';
set(lgnd,'Interpreter','latex')


