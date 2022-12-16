clc; clear; close all;

%% assegno i dati
LL = 4; %m
DD = .01; %m
II = 1e3; %A
kk = 350; %W/mK
hh = 500; %W/m^2K
roel = 1.75e-8; %hom*m
Tb = 77; %K

SE = pi*DD*LL; % sup ext
SI = pi*DD^2/4; % cross section
VV = SI*LL; % volume

RR = roel*LL/(pi*DD^2/4); %hom
qqq = (RR*II^2)/VV + SE*hh/VV*Tb; %W/m^3

%% creazione griglia
dx = 5e-3;
xx = (0:dx:LL)';
nn = length(xx);

%% assembraggio matrice AA
s_diag = ones(nn,1);
c_diag = (-2-SE*hh/kk/VV*dx^2)*ones(nn,1);
i_diag = ones(nn,1);

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% assembraggio vettore bb
bb = -qqq*dx^2/kk*ones(nn,1);

%% coco
AA(1,1) = 1; AA(1,2) = 0; bb(1) = Tb; % dirichlet
AA(end,end-1) = 0; AA(end, end) = 1; bb(end) = Tb; % dirichlet

%% soluzione numerica
TT = AA\bb;

%% plotto
plot(xx, TT, 'linewidth', 2)
legend('soluzione numerica')
grid on
xlabel('x [m]')
ylabel('T [K]')
