clc; 
clear; 
close all;

%% dati
% fisici
Tco2 = 300;         % °C
hh = 50;            % W/m^2K
qvol = 8e4;         % W/m^3
kbarra = 350;       % W/mK
kss = 50;           % W/mK


% topologia
d1 = .2;            % m
d2 = .23;           % m
LL = .5;            % m

surfi1 = pi*d1^2/4;             % m^2
surfi2 = pi*(d2^2-d1^2)/4;      % m^2
surfe = pi*(d1+d2)*LL/2;          % m^2 COSTANTE

vol1 = surfi1*LL;           % m^3
vol2 = surfi2*LL/2;         % m^3
voltot = vol1+vol2;         % m^3 COSTANTE

%% griglia
dx = 1e-4;                  % dalla considerazione dell'andamento errori !     
x1 = (0:dx:LL/2)';    
x2 = (LL/2+dx:dx:LL)';
xx = [x1;x2];
n1 = length(x1);
n2 = length(x2);
nn = length(xx);

keq = (kbarra*(surfi1)+kss*(surfi2))/(surfi1+surfi2);
kk=kbarra*(xx<=LL/2)+keq*(xx>LL/2);

%% AA & bb
diag_i = ones(nn,1);
diag_p = -2+surfe/voltot*dx^2./-kk*hh;
diag_s = ones(nn,1);

AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);

bb = qvol*dx^2./-kk+surfe/voltot*hh*Tco2*dx^2./-kk;

%% coco & TT
% robin
AA(1,1) = -1+hh*dx/-kbarra;
AA(1,2) = 1;
bb(1) = hh/-kbarra*Tco2*dx;

% continuità
AA(n1,n1-1) = -kbarra;
AA(n1,n1) = kbarra+keq;
AA(n1,n1+1) = -keq;
bb(n1) = 0;

% robin
AA(end,end-1) = -1;
AA(end,end) = 1-hh/-keq*dx;
bb(end) = hh/-keq*-Tco2*dx;

TT = AA\bb;

%% plotto
figure(1)
plot(xx*1e2,TT,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Lunghezza (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)