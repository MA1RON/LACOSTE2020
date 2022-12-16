clc; clear; close all;

%% dati
ee = 35e-3; %m
ss = 30e-3; %m
dd = 25e-3; %m
rin = dd/2;
kk = 17; %W/mK
Tw = 450; %K
hhin = 1e3; %W/m^2K
Ta = 298.15; %K
hhout = 1e2; %W/m^2K

surf=3*sqrt(3)*(ee/2)^2/2; % Area
DD=sqrt(4*surf/pi); % D idraulico
rout = DD/2;

%% creazione griglia
dr = 1e-4;
rr = (rin:dr:rout)';
nn = length(rr);

%% assembraggio matrice AA
sub_diag = 1-dr/2./rr;
principale = -2*ones(nn,1);
super_diag = 1+dr/2./(rr);
BB = [ [sub_diag(2:end);0] principale [0;super_diag(1:end-1)] ];

AA = spdiags(BB, -1:1, nn, nn);

%% assembraggio vettore bb
bb = zeros(nn,1);

%% coco
% robin
AA(1,1) = dr/kk*hhin+1; 
AA(1,2) = -1; 
bb(1) = dr*hhin/kk*Tw;
% robin
AA(end,end-1) = -1; 
AA(end, end) = dr/kk*hhout+1; 
bb(end) = dr*hhout/kk*Ta; 

%% soluzione numerica
TT = AA\bb;

%% soluzione analitica || C1 C2 T1 T2
%metodo matriciale
AAan = [-log(rin), -1, 1, 0;
        -log(rout), -1, 0, 1;
        1/rin, 0, hhin/-kk, 0;
        1/rout, 0, 0, hhout/kk];
bban = [0, 0, hhin*Tw/-kk, hhout*Ta/kk]';
TNan = AAan\bban;
C1 = TNan(1); C2 = TNan(2); T1an = TNan(3); T2an = TNan(4);
Tan= @(r) C1*log(r)+C2;
%}

%a manina santa
C1=hhout*(Tw-Ta)/(-kk/rout-hhout*log(rout)-hhout/hhin*kk/rin+hhout*log(rin));
C2=Tw+(kk/rin-hhin*log(rin))*C1/hhin;
%}
Tan_by_prof= @(r) C1*log(r)+C2;

%% plotto
plot(rr, TT, 'linewidth', 2)
hold on
plot(rr, Tan_by_prof(rr),'r-.', 'linewidth', 2)
hold on
plot(rr, Tan(rr),'k-.', 'linewidth', 2)
legend('soluzione numerica', 'analitica a manina', 'analitica matriciale')
set(gca,'fontsize',22)
grid on
xlim([rin rout])
xlabel('r [m]')
ylabel('T [K]')