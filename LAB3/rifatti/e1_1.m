clc; clear; close all;

%% dati
delta = 50e-3; %m
kk = 5; %W/mK
qqq = 5e5; %W/m^3
T0 = 300; %K

%% griglia
dx = 5e-3; %m
xx = (0:dx:delta)';
nn = length(xx);

%% AA
s_diag = ones(nn,1);
c_diag = -2*ones(nn,1);
i_diag = ones(nn,1);

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% bb
bb = ones(nn,1)*qqq/-kk*dx^2;

%% coco
AA(1,1) = -1; AA(1,2) = 1; b(1) = 0; %neumann
AA(end, end-1) = 0; AA(end, end) = 1; bb(end) = T0; %dirichlet

%% TT
TT = AA\bb;

%% analitica
Tan = @(x) qqq/kk/2*(delta^2-x.^2)+T0;

%% show
plot(xx, TT, 'r-', 'linewidth', 2);
hold on
plot(xx, Tan(xx), 'b-.', 'linewidth', 2);
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('x [m]')
ylabel('T [K]')