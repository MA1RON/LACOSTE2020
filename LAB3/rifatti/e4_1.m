clc; clear; close all;

%% assegno i dati
delta = 50e-3; %m
kk = 5; %W/mK
qvol = 5e5; %W/m^3
T0 = 300; %K
Tf = 273; %K
hh = 100; %W/m^2K

%% creazione griglia
dx = 5e-5;
xx = (0:dx:delta)';
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
AA(1,1) = -1-2*hh*dx/kk; AA(1,2) = 1; bb(1) = -hh*Tf*dx*2/kk; % robin 
AA(end,end-1) = 0; AA(end, end) = 1; bb(end) = T0; % dirichlet

%% soluzione numerica
TT = AA\bb;

%% soluzione analitica
C2 = (T0+qvol*delta^2/2/kk+hh*Tf/kk*delta)/(1+hh*Tf/kk*delta);
C1 = hh*(Tf-C2)/-kk;
Tan = @(x) -qvol/kk*x.^2+C1*x+C2;

%% plotto
plot(xx, TT, 'linewidth', 2)
hold on
plot(xx, Tan(xx), 'r-.', 'linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('x [m]')
ylabel('T [K]')
