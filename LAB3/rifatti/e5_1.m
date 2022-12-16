clc; clear; close all;

%% dati
Tw = 400; %K
ri = 10e-3; %m
ro = 12e-3; %m
kk = 0.035; %W/mK
hi = 100; %W/m^2K
Ta = 300; %K
ho = 10; %W/m^2K

%% griglia
dr = 5e-5; %m
rr = (ri:dr:ro)';
nn = length(rr);

%% AA
s_diag = -dr + 2*rr;
c_diag = -4*rr;
i_diag = dr + 2*rr;

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% bb
bb = zeros(nn,1);

%% coco
%AA(1,1) = 1; AA(1,2) = 0; bb(1) = 400; % dirichlet
AA(1,1) = -1 + dr*hi/-kk; AA(1,2) = 1; bb(1) = dr*Tw*hi/-kk; % robin
%AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = 300; % dirichlet
AA(end, end-1) = -1; AA(end,end) = 1 - dr*ho/-kk; bb(end) = -dr*Ta*ho/-kk; % robin

%% TT
TT = AA\bb;

%% analitica || C1 C2 T1 T2
C1 = (-Ta+Tw)/kk/(1/ri/hi-1/ro/ho+log(ro/ri));
C2 = Tw-C1*log(ri)+kk/ri/hi*C1;
TTan = @(r) C1*log(r) + C2;
  
%% show
plot(rr, TT, 'r-', 'linewidth', 2)
hold on
plot(rr, TTan(rr), 'b-.', 'linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('r [m]')
ylabel('T [K]')