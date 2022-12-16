clc; clear; close all;
%% dati
% dominio
% x -> Tinf | 0 - barra (ba) - contatto - barra&guarnizione (ss) - end | Tinf
% geometrici
DDint = .2;
LL = .5;
DDext = .23;

surfint = DDint^2*pi/4;
volint = surfint*LL;

surfext = (DDext^2-DDint^2)*pi/4;
volext = surfext*LL/2;

voltot = volint + volext;
Awet = pi*DDint*LL/2 + pi*DDext*LL/2 + surfint + (surfext+surfint) + surfext;

% fisici
Tinf = 300;
hh = 50;
qint = 80e3;
kkba = 350;
kkss = 50;

kkave = (kkba*surfint + kkss*surfext)/(surfint+surfext);

%% griglia spaziale
dx = 1e-3;
xx = (0:dx:LL)';
nn = length(xx);
nncont = round(nn/2); % contatto

% dati dipendenti da x
kk = kkba*(xx <= xx(nncont)) + kkave*(xx > xx(nncont));

%% sistema
diag_si = ones(nn,1);
diag_p = -2 - hh./kk*dx^2*Awet/voltot;
AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
bb = -(qint + hh*Tinf*Awet/voltot)./kk*dx^2;

% bb = -qint*dx^2./kk - hh*Tco2*As/VV./kk*dx^2;

% --- condizioni al contorno ---
% robin iniziale
AA(1,1) = 1+hh/kkba*dx;
AA(1,2) = -1;
bb(1) = hh/kkba*dx*Tinf;

% contatto
AA(nncont,nncont-1) = kkba/kkave;
AA(nncont,nncont) = -1-kkba/kkave;
AA(nncont,nncont+1) = 1;
bb(nncont) = 0;

% robin finale
AA(end,end-1) = -1;
AA(end,end) = 1+hh/kkave*dx;
bb(end) = hh/kkave*dx*Tinf;

% --- risolvo il sistema ---
TT = AA\bb;

%% post production
figure
plot(xx*1e2,TT,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Lunghezza (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)