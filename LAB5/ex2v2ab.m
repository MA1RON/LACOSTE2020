clc; clear; close all;
%% dati
% dominio
% Tint | 0 - Cu - contatto - is - end | Text
% geometrici
DDint = 30e-3;
DDcont = 50e-3;
DDext = 52e-3;
LL = 1; % come se non ci fosse, mi aiuta solo il ragionamento

surfcu = pi*(DDcont^2-DDint^2)/4;
surfis = pi*(DDext^2-DDcont^2)/4;
surftot = pi*(DDext^2-DDint^2)/4;
volcu = surfcu*LL;
vol = surftot*LL;

Awetint = pi*DDint*LL;
Awetext = pi*DDext*LL;

% fisici
kcu = 150;
roelcu = 2e-8;
II = 20e3;
Tint = 77;
hhint = 1200;
kis = .5;
Text = 25+273.15;
hhext = 15;

RR = roelcu/surfcu*LL;
qqq = RR*II^2/volcu;

cpcu = 350;
rocu = 8900;
cpis = 1800;
rois = 2500;

%% costo computazionale
drv = [1e-6 2e-6 3e-6 5e-6 1e-5 2e-5 3e-5 5e-5 1e-4 2e-4 3e-4 5e-4];
NN = zeros(size(drv));
costo = zeros(size(drv));
for jjerr = 1:length(drv)
    %% griglia spaziale
    dr = drv(jjerr);
    rr = [(DDint/2:dr:DDcont/2)';(DDcont/2+dr:dr:DDext/2)'];
    nncont = length(DDint/2:dr:DDcont/2);
    nn = length(rr);

    kk = kcu*(rr <= rr(nncont)) + kis*(rr > rr(nncont));

    %% creo il sistema
    diag_s = .5+rr/dr;
    diag_p = -2*rr/dr;
    diag_i = -.5+rr/dr;
    AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1, nn,nn);
    bb = qqq/-kcu*dr*rr.*(rr<=rr(nncont));

    % --- condizioni al contorno ---
    AA(1,2) = -1;
    AA(1,1) = 1+hhint/kk(1)*dr;
    bb(1) = hhint/kk(1)*dr*Tint;

    AA(nncont,nncont-1) = -1;
    AA(nncont,nncont) = 1+kis/kcu;
    AA(nncont,nncont+1) = -kis/kcu;
    bb(nncont) = 0;

    AA(end,end-1) = 1;
    AA(end,end) = -1-hhext/kk(end)*dr;
    bb(end) = -hhext*Text/kk(end)*dr;

    %% risolvo il sistema
    tic
    TT = AA\bb;
    costo(jjerr) = toc;
    NN(jjerr) = nnz(AA);
    
    if jjerr == 1 % soluzione più precisa
        TTref = TT;
        rrref = rr;
        nncontref = nncont;
    end
end

%% post prodution
% andamento temperatura
figure
subplot(1,2,1)
plot(rrref*1e3,TTref,'b-.','linewidth',2)
hold on
plot([rrref(nncontref) rrref(nncontref)]*1e3,[min(TTref) max(TTref)],'k--','linewidth',2)
legend('Temperatura','Superficie di contatto','location','northwest')
grid on
box on
xlabel('Raggio [mm]')
ylabel('Temperatura [K]')
title('Stazionario 1D')
set(gca,'fontsize',16)

subplot(1,2,2)
loglog(NN,costo,'b-.*','linewidth',2)
grid on
box on
xlabel('Non zero elements [-]')
ylabel('Costo computazionale [s]')
title('Costo stazionario 1D')
set(gca,'fontsize',16)

%% conservazione dell'energia
qvol = roelcu*II^2*LL/surfcu/volcu; % W/m^3
Eint = qvol*volcu % W
Eout = abs(hhext*(Text-TT(end)))*Awetext...
       + abs(hhint*(Tint-TT(1)))*Awetint % W