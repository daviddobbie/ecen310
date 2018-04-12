%% David Dobbie
% 300340161
% Assignment 2, Question 4, ECEN 310

close all;
clear;
clf;
clc;


% creates the noise for an exact simulated Pm

N0 = 1;
Eb = 1;

count = 10000;

% init const. points on 2D plane
bit0 = [sqrt(Eb)];
bit1 = [-sqrt(Eb)];
bit2 = [sqrt(Eb)*3];
bit3 = [-sqrt(Eb)*3];

decisionPoints = [bit0; bit1; bit2; bit3];  

SNR = Eb/N0; %since all equal, can do this

% init 0 mean noise on 2D plane
noise = normrnd(0, N0/2, [count 1]);

data = [bit0.*(ones(count/4,2)); bit1.*(ones(count/4,2)); 
    bit2.*(ones(count/4,2)); bit3.*(ones(count/4,2))] + noise;

figure(1)
hold on
ax = gca;
ax.YAxisLocation = 'origin';  % setting y axis location to origin
ax.XAxisLocation = 'origin';  % setting y axis location to origin
grid on
histogram(data,200)

xlabel('Signal Energy in terms of $\sqrt E_b$','interpreter','latex');
ylabel('Number of bits arrived','interpreter','latex');


% decision making finction, for exact simulated Pm


% plots the analytical expression found in class


% plots the union bound accounding to the min diatance approx


% plots the nbound using min dist approx and bound on q func