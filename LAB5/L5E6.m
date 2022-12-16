clc; 
clear; 
close all;

%% dati
% fisici
Tin = 25;       % °C
hin = 2e3;      % W/m^2K
Tout = 15;      % °C
hout = 5;       % W/m^2K

Inom = 16e3;    % A

kss = 1.5;      % W/mK
kcu = 30;       % W/mK
kis = .21;      % W/mK

ross = 5e-6;    % hom*m
rocu = 1.8e-8;  % hom*m
rois = 1e-2;    % hom*m

% topologia
r1 = 10e-3;     % m
r2 = 14e-3;     % m
r3 = 17e-3;     % m
r4 = 18.5e-3;   % m

surfiss = pi*(r2^2-r1^2);                           % m^2
surficu = pi*(r3^2-r2^2);                           % m^2
surfiis = pi*(r4^2-r3^2);                           % m^2
surfitot = surfiss + surficu + surfiis;             % m^2

% partitore di corrente
Rss = ross/surfiss;                                 % hom/m
Rcu = rocu/surficu;                                 % hom/m
Ris = rois/surfiis;                                 % hom/m

Iss = Inom*(Rcu*Ris)/(Rss*Rcu+Rcu*Ris+Ris*Rss);     % A
Icu = Inom*(Rss*Ris)/(Rss*Rcu+Rcu*Ris+Ris*Rss);     % A
Iis = Inom*(Rcu*Rss)/(Rss*Rcu+Rcu*Ris+Ris*Rss);     % A

%% griglia
dr = 1e-4;
x1 = (r1:dr:r2)';
n1 = length(x1);
x2 = (r2+dr:dr:r3)';
n2 = length(x2);
x3 = (r3+dr:dr:r4)';
n3 = length(x3);
rr = [x1;x2;x3];
nn = length(rr);

% parametri vettoriali
II = Iss*(rr<=r2) + Icu*(rr>r2 & rr<=r3) + Iis*(rr>r3);
kk = kss*(rr<=r2) + kcu*(rr>r2 & rr<=r3) + kis*(rr>r3);
RR = Rss*(rr<=r2) + Rcu*(rr>r2 & rr<=r3) + Ris*(rr>r3);
ro = ross*(rr<=r2) + rocu*(rr>r2 & rr<=r3) + rois*(rr>r3);
surfi = surfiss*(rr<=r2) + surficu*(rr>r2 & rr<=r3) + surfiis*(rr>r3);

qqq = ro.*II.^2./surfi.^2;

%% AA & bb
diag_s = dr+2*rr;
diag_p = -4*rr;
diag_i = -dr+2*rr;
AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);

bb = qqq./-kk*2*dr^2.*rr;

%% coco & TT
% robin
AA(1,1) = -1+hin/-kss*dr;
AA(1,2) = 1;
bb(1) = hin/-kss*dr*Tin;

% continuità
AA(n1,n1-1) = kss;
AA(n1,n1) = -kss-kcu;
AA(n1,n1+1) = kcu;
bb(n1) = 0;

% continuità
n2 = n2 + n1;           % OCCHIO
AA(n2,n2-1) = kcu;
AA(n2,n2) = -kcu-kis;
AA(n2,n2+1) = kis;
bb(n2) = 0;

% robin
AA(end,end-1) = -1;
AA(end,end) = 1+hout/kis*dr;
bb(end) = hout/kis*dr*-Tout;

TT = AA\bb;

%% plotto
figure(1)
plot(rr*1e3,TT,'linewidth',2)
grid on
box on
title('Distribuzione di temperatura')
xlabel('Raggio (mm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)