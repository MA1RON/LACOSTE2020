clc; clear; close all;

%% dati
Ta = 5; %°C
hh = 100; %W/m^2K
qq = 50; %W/m^2
rin = .1; %m
rout = .15; %m
kk = .2; %W/mK

%% griglia
dr = 1e-3;
rr = (rin:dr:rout)';
nn = length(rr);

%% AA
diag_s=1-dr./(rr);
diag_p=-2*ones(nn,1);
diag_i=1+dr./(rr);
BB = [[diag_s(2:end);0] diag_p [0;diag_i(1:end-1)]];
AA = spdiags(BB, -1:1, nn, nn);

%% bb
bb = zeros(nn,1);

%% coco
%neumann non omo
AA(1,2) = -1;
AA(1,1) = 1;
bb(1) = qq/kk*dr;

%robin
AA(end,end-1) = -1;
AA(end,end) = 1 + hh*dr/kk;
bb(end) = hh*Ta*dr/kk;

%% TT
TT = AA\bb;

%% analitica
C1 = -qq*rin^2/kk;
T2 = -kk/hh*C1/rout^2+Ta;
C2 = T2+C1/rout;
Tan = @(r) -C1./r + C2;

%% plotto
plot(rr, TT, 'linewidth', 2)
hold on
plot(rr, Tan(rr),'r-.', 'linewidth', 2)
legend('soluzione numerica', 'analitica')
grid on
set(gca, 'fontsize', 18)
xlabel('Raggio [m]')
ylabel('Temperatura [°C]')