
clear all
clc 



set(0,'defaulttextInterpreter','latex')

EGC = 1 + pi/4;
SC = 3/2;
MRC = 2;

b= 0: 0.1 : 5;


meanSNR = [1 2.5 5];








figure(1)
clf

hold on
SSC = 1 + exp(-b).*( meanSNR(1) + b  -1);
plot(b,SSC)
plot([0 max(b)], [EGC EGC]);
plot([0 max(b)], [SC SC]);
plot([0 max(b)], [MRC MRC]);
hold off


xlabel('$\frac{\gamma_T}{\bar{\gamma}}$ Threshold-Mean Ratio')
ylabel('Combined Gain $\frac{\bar{\gamma_{C}}}{\bar{\gamma}}$')

lgf = legend('SSC $\bar{\gamma} =1$', 'EGC','SC','MRC');
set(lgf,'Interpreter','latex');

ylim([1 2.5])