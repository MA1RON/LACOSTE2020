clc; clear; close all; 

%% dati
kk = .15; %W/mK
ri = 10e-3; %m
ro = 25e-3; %m
Tw = 15; %°C
hh = 1000; %W/m^2K
qq = -2e3; %W/m^2 OCCHIO AL SEGNO

%% griglia
dr = 5e-4; %m
rr = (ri:dr:ro)';
nn = length(rr);

%% AA
s_diag = dr + 2*rr;
c_diag = -4*rr;
i_diag = -dr + 2*rr;

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% bb
bb = zeros(nn, 1);

%% coco
AA(1,1) = -1-dr*hh/kk; AA(1,2) = 1; bb(1) = -dr*hh*Tw/kk; % robin
AA(end, end-1) = -1; AA(end, end) = 1; bb(end) = -qq*dr/kk; % neumann

%% TT
TT = AA\bb;

%% analitica
A = -ro*qq/kk;
T1 = Tw+A/ri*kk/hh;
B = T1-A*log(ri);
Tan = @(r) A*log(r)+B;

%% show
plot(rr, TT, 'r-', 'linewidth', 2)
hold on
plot(rr, Tan(rr), 'b-.', 'linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('r [m]')
ylabel('T [°C]')