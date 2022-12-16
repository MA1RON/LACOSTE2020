clc; clear; close all;

%% assegno i dati
LL = 50e-3; %m
kk = 5; %W/mK
qvol = 5e5; %W/m^3
T0 = 300; %K

%% creazione griglia
dx = 5e-3;
xx = (0:dx:LL)';
nn = length(xx);

%% assembraggio matrice AA
sub_diag = ones(nn,1);
principale = -2*ones(nn,1);
super_diag = ones(nn,1);
BB = [sub_diag principale super_diag];

AA = spdiags(BB, -1:1, nn, nn);

%% assembraggio vettore bb
bb = -qvol*ones(nn,1)*dx^2/kk;

%% coco
AA(1,1) = -1; AA(1,2) = 1; bb(1) = 0; % neumann omo
AA(end,end-1) = 0; AA(end, end) = 1; bb(end) = T0; % dirichlet

%% soluzione numerica
TT = AA\bb;

%% soluzione analitica
Tan = @(x) -qvol/kk/2*(x.^2-LL^2) + T0;

%% plotto
plot(xx, TT, xx, Tan(xx), 'linewidth', 3)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('x [m]')
ylabel('T [K]')
set(gca,'fontsize',18)
