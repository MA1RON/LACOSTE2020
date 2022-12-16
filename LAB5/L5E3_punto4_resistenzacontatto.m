clc; clear; close all;

%% dati
tu = 2e-3; %m
td = 1e-3; %m
Tu = 20; %°C
Td = 70; %°C
hh = 15; %W/m^2K
ku = 1; %W/mK
kd = 15; %W/mK

tis = 0.01e-3; %5e-4 da testo m
Ris = 2.5e-4; %2.5e-5 da testo Km^2/W
kguarn = tis/Ris; %W/mK

%% griglia {(up) 0->robin->i-> is ->s->dirichlet->end (down)}
dx = 1e-6; %m
xu = (0:dx:tu)';
nu = length(xu);
xis = (tu+dx:dx:tu+tis)';
nis = length(xis);
xd = (tu+tis+dx:dx:tu+td+tis)';
xx = [xu;xis;xd];
nn = length(xx);

%% AA & bb
my_diag = ones(nn,1);
AA = spdiags([my_diag -2*my_diag my_diag], -1:1, nn, nn);
bb = zeros(nn,1);

%% coco & TT
%robin
AA(1,1) = hh/ku*dx+1;
AA(1,2) = -1;
bb(1) = hh*Tu*dx/ku;

%continuità
AA(nu,nu-1) = -ku;
AA(nu,nu) = ku + kguarn;
AA(nu,nu+1) = -kguarn;
bb(nu) = 0;

%continuità
AA(nu+nis,nu+nis-1) = -kguarn;
AA(nu+nis,nu+nis) = kd + kguarn;
AA(nu+nis,nu+nis+1) = -kd;
bb(nu+nis) = 0;

%dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Td;

TT = AA\bb;

%% plotto
plot(xx,TT, 'b-', 'linewidth', 2)
hold on
plot([2e-3 2e-3], [70 68], 'k-.', 'linewidth', 2)
legend('Temperatura', 'Contatto', 'location', 'south')
title('Distribuzione di temperatura')
xlabel('Spessore (mm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)
box on
grid on
xlim([xx(1) xx(end)])