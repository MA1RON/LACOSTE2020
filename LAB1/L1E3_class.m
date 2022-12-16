clear; close all; clc;

AA = [2:5;0 2 6 0;0 0 2 9;.5:.5:2];

%% 3.1
BB = [max(AA); min(AA')];

%% 3.2
find(BB)

%% 3.3
AA*BB'

%% 3.4
AA\BB'

%% 3.5
AA^2
AA.^2