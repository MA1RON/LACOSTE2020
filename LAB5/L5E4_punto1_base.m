clc; clear; close all;

%% dati
LL = .5; %m
bb = .1; %m
aa = .02; %m
kk = 350; %W/mK
ro = 8900; %kg/m^3
cp = 350; %J/kgK
THe = 4.5; %K
hh = 500; %W/m^2K
II = 7e3; %A
roel0 = 1.75e-8; %hom*m
L0 = .56; %m
roel = @(z) roel0*cos(pi*z/L0).^2; %hom*m

% parametri topologici
surfi = aa*bb; %m^2
surfe = 2*(aa+bb)*LL; %m^2
vol = surfi*LL; %m^3

%% griglia con origine a metà
dz = 1e-4;
zz = (0:dz:LL/2)';
nn = length(zz);

%% qqq = qqel + qqconv(non f(T))
qqel = roel(zz)*LL/surfi*II^2/vol;
qqconv = hh*THe*surfe/vol;
qqq = qqel + qqconv;

%% AA & bb
diag_s = ones(nn,1);
diag_p = -2*ones(nn,1)+hh*dz^2*surfe/vol/-kk;
diag_i = ones(nn,1);
AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);

bb = qqq*dz^2/-kk;

%% coco & TT
% neuman omo (simmetria)
AA(1,1) = 1;
AA(1,2) = -1;
bb(1) = 0;

% dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = THe;

TT = AA\bb;

%% plotto
plot(zz,TT, 'b-', 'linewidth', 2)
title('Distribuzione di temperatura')
grid on
xlabel('z [m]')
ylabel('T [K]')
set(gca, 'fontsize', 18)
