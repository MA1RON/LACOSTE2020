clc; clear; close all; 

%% dati
delta = 50e-3; %m
kk = 5; %W/mK
qqq = 5e5; %W/m^3
T0 = 300; %K
Tf = 273; %K
hh = 100; %W/m^2K

%% griglia
dx = 5e-3; %m
xx = 0:dx:delta;
nn = length(xx);

%% AA
s_diag = ones(nn,1);
c_diag = -2*ones(nn,1);
i_diag = ones(nn,1);

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% bb
bb = qqq/-kk*ones(nn,1)*dx^2;

%% coco
AA(1,1) = -1+hh/-kk*dx; AA(1,2) = 1; bb(1) = hh/-kk*Tf*dx; % robin
AA(end, end-1) = 0; AA(end,end) = 1; bb(end) = T0; % dirichlet

%% TT
TT = AA\bb;

%% analitica
C1 = ((T0-Tf)/delta+qqq*delta/2/kk)/(1+kk/hh/delta);
C2 = kk*C1/hh+Tf;
Tan = @(x) qqq/-kk*x.^2/2 +C1*x+C2;

plot(xx,TT, 'r-', 'linewidth', 2)
hold on
plot(xx,Tan(xx), 'b-.', 'linewidth', 2)
grid on
legend('soluzione numerica', 'soluzione analitica')
xlabel('x [m]')
ylabel('T [k]') 