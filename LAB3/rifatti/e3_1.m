clc; clear; close all; 


%% dati
DD = .01; %m
LL = 4; %m
II = 1e3; %A
Tb = 77; %K
hh = 500; %W/m^2K
roel = 1.75e-8; %hom*m
kk = 350; %W/mK

%% griglia
dx = 5e-3; %m
xx = 0:dx:LL;
nn = length(xx);

%% qqq
SE = pi*DD*LL; %m^2
SI = pi*DD^2/4; %m^2
VV = SI*LL; %m^3
RR = roel*LL/SI; %hom
qqq = RR*II^2/VV + hh*SE*Tb/VV; %W/m^3

%% AA
s_diag = ones(nn,1);
c_diag = ones(nn,1)*(-2-dx^2*hh/kk*SE/VV);
i_diag = ones(nn,1);

AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);

%% bb
bb = (-qqq*dx^2/kk)*ones(nn,1);

%% coco
AA(1,1) = 1; AA(1,2) = 0; bb(1) = Tb;
AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = Tb;

%% TT
TT = AA\bb;

%% show
plot(xx, TT, 'r-', 'linewidth', 2)
grid on
legend('soluzione numerica', 'location', 'south')
xlabel('x [m]')
ylabel('T [K]')