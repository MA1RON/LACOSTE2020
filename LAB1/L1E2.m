clear; close all; clc;

AA = 2:2:8;
BB = logspace(2,11,4);

%% 2.1
CC = log10(BB);

%% 2.2
DD = [AA; CC];

rr = size(DD,1);
cc = size(DD,2);

%% 2.3
DD = [DD; ones(1,cc) * 3; zeros(1,cc)];

%% 2.4
EE = DD * DD';
RR = sum(EE,1);
KK = sum(EE,2);
TT = RR~=0 & KK'>100;