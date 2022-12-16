clc; clear; close all;
%% dati
% topologia
LL = 4;                     % m
DD = 33.7e-3;               % m
tt = 2.9e-3;                % m

rou = DD/2;                 % m
rin = rou-tt;               % m

% fisici
kac = 13;                   % W/mK
Tw = 100+273.15;            % K
Pr = 0.7;                   % -
Re = 1e5;                   % -
khe = @(p, T) 2.682e-3*(1 + 1.123e-3*p).*T.^(0.71*(1-2e-4*p));  % W/mK
Nu_function = @(Pr, Re) 0.023*Re^0.8*Pr^0.7;                    % Dittus-Boelter alla T media !
Nu = Nu_function(Pr, Re);
Thein = 20+273.15;          % K
Theou = 100+273.15;         % K
Theave = (Thein+Theou)/2;   % K
phe = 5;                    % bar
hhe = khe(phe,Theave)*Nu/2/rou; % W/m^2K

%% griglia radiale
dr = 1e-4;                  % m
rr = (rin:dr:rou)';
nn = length(rr);

%% AA & bb
diag_i = -dr+2*rr;
diag_p = -4*rr;
diag_s = dr+2*rr;
AA = spdiags([[diag_i(2:end);0], diag_p [0;diag_s(1:end-1)]], -1:1, nn, nn);

bb = zeros(nn,1);

%% coco & TT
% robin int
AA(1,1) = -1+hhe/-kac*dr;
AA(1,2) = 1;
bb(1) = hhe*dr/-kac*Thein; % prendo l'entrata al tubo con sbalzo maggiore

% dirichlet ext
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tw;

TT = AA\bb;

%% plotto
plot(rr*1e3,TT-273.15, 'b-', 'linewidth', 2)
grid on
set(gca, 'fontsize', 18)
box on
xlabel('Raggio [mm]')
ylabel('Temperatura [°C]')
title('Andamento radiale')
