clc; clear; close all;

%% assegno i dati
din = 20e-3; %m
dout = 50e-3; %m
kk = .15; %W/mK
qsurf = 2e3; %W/m^3
Tw = 15; %°C
hh = 1e3; %W/m^2/K

%% creazione griglia
dr = 1e-3;
rr = (din/2:dr:dout/2)';
nn = length(rr);

%% assembraggio matrice AA
s_diag = 1-dr/2./rr;
p_diag = -2*ones(nn,1);
i_diag = 1+dr/2./(rr);
BB = [ [s_diag(2:end);0] p_diag [0;i_diag(1:end-1)] ];

AA = spdiags(BB, -1:1, nn, nn);

%% assembraggio vettore bb
bb = zeros(nn,1);

%% coco
% robin conv
AA(1,1) = dr/kk*hh+1; 
AA(1,2) = -1; 
bb(1) = dr*hh/kk*Tw; 
% neumann non omo
AA(end,end-1) = 1; 
AA(end, end) = -1; 
bb(end) = -dr/kk*qsurf; 

%% soluzione numerica
TT = AA\bb;

%% soluzione analitica
C1 = dout/2/-kk*qsurf; 
T1 = kk*C1/din*2/hh+Tw;
C2 = T1-C1*log(din/2);
%Tan = @(r) Tw + qsurf*(dout/2)*(log((din/2)/r)/kk+1/(hh*(din/2)));
Tan= @(r) qsurf*dout/2/kk*(log(r/(din/2))+kk/(din/2*hh))+Tw;

%% plotto
plot(rr, TT, 'linewidth', 2)
hold on
plot(rr, Tan(rr),'r-.', 'linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('r [m]')
ylabel('T [°C]')
