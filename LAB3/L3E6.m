clc; clear; close all;

%% assegno i dati
tb = .25; %m
ti = .1; %m
kb = .55; %W/mK
ki = .04; %W/mK
Tin = 20; %°C
Text = -5; %°C
hext = 25; %W/m^2K

tt = tb+ti;
kk = (kb*tb+ki*ti)/tt;

%% creazione griglia
dx = 5e-5;
xx = (0:dx:tt)';
nn = length(xx);

%% assembraggio matrice AA
diag_i = ones(nn,1);
diag_p = -2*ones(nn,1);
diap_s = ones(nn,1);
AA = spdiags([diag_i diag_p diap_s], -1:1, nn, nn);

%% assembraggio vettori bb
bb = zeros(nn, 1);

%% coco
% robin
AA(1,1) = -1+hext*dx/-kk; 
AA(1,2) = 1; 
bb(1) = hext*Text*dx/-kk;
% dirichlet
AA(end,end-1) = 0; 
AA(end, end) = 1; 
bb(end) = Tin;

%% soluzione numerica
TT = AA\bb;

%% plotto
plot(xx, TT, 'linewidth', 2)
legend('soluzione numerica')
grid on
xlabel('x [m]')
ylabel('T [K]')
