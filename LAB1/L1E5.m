clear; close all; clc;

QQ = diag([3 4], -2) + diag(pi/2:pi/2:3*pi/2, -1) + diag(-12:5:3) + diag(3*pi/2:pi/2:5*pi/2, 1) + diag([-2 -1], 2);
QQ(end) = 1;

%% 5.1
QQ(:,end) = QQ(:,1);

%% 5.2
QQ(:,3) = QQ(2,:);

%% 5.3
NN = [QQ(1,:)' QQ(:,2)];

%% 5.4
QQ'*NN