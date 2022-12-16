clc; 
%clear; 
close all;

%% dati
% fisici
Tco2 = 300;         % °C
hh = 50;            % W/m^2K
qvol = 8e4;         % W/m^3
kbarra = 350;       % W/mK
kss = 50;           % W/mK

% topologia
d1 = .2;            % m
d2 = .23;           % m
LL = .5;            % m

surfi1 = pi*d1^2/4;             % m^2
surfi2 = pi*(d2^2-d1^2)/4;      % m^2
surfe = pi*(d1^2+d2^2)/4;       % m^2 COSTANTE (perché entrambi??)

vol1 = surfi1*LL;           % m^3
vol2 = surfi2*LL/2;         % m^3
voltot = vol1+vol2;         % m^3 COSTANTE

%% griglia
dr = 1e-4;
rr = (0:dr:d1/2)';
nn = length(rr);

kk = kbarra;

%% AA & bb
diag_s = dr+2*rr;
diag_p = -4*rr + 2*rr*surfe/voltot*hh/-kk*dr^2;
diag_i = -dr+2*rr;
AA = spdiags([[diag_i(2:end); 0] diag_p [0; diag_s(1:end-1)]], -1:1, nn, nn);

bb = 2*dr^2*rr*(qvol+hh*Tco2*surfe/voltot)/-kk;

%% coco & TT
% neuman omo centrale
AA(1,1) = 1;
AA(1,2) = -1;
bb(1) = 0;

% robin esterna
AA(end,end-1) = -1;
AA(end,end) = 1 + hh*dr/kk;
bb(end) = hh*dr/kk*Tco2;

TT = AA\bb;

%% plotto
figure(1)
plot(rr*1e2,TT,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Raggio (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)