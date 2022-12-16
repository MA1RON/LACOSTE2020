clc; clear; close all;
%% dati
% geometrici
DDint = 1e-2;
DDext = 2e-2;

% fisici
kk = .1;
qqq = 2e6;
Text = 300;
hhext = 10;
%% verifico errore
drv = [1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3];
err = zeros(size(drv));
errcond = zeros(size(drv));
Tmaxv = zeros(size(drv));
for jjerr = 1:length(drv)
    %% griglia
    dr = drv(jjerr);
    rr = (DDint/2:dr:DDext/2)';
    nn = length(rr);

    %% costruisco il sistema - adimensionale!
    diag_s = .5+rr/dr;
    diag_p = -2*rr/dr;
    diag_i = -.5+rr/dr;
    AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1,nn,nn);
    bb = -qqq/kk*rr*dr;

    % --- condizioni al contorno ---
    % interno adiabatico
    AA(1,1) = 1;
    AA(1,2) = -1;
    bb(1) = 0;

    %{
    % robin esterno - polinomiale
    AA(end,end-2) = 3;
    AA(end,end-1) = -4;
    AA(end,end) = 1-hhext*2*dr/kk;
    bb(end) = hhext*2*dr*-Text/kk;
    %}
    % esterno convettivo
    AA(end,end-1) = 1;
    AA(end,end) = -1-hhext*dr/kk;
    bb(end) = hhext*dr*-Text/kk;
    %}

    %% risolvo il sistema
    TT = AA\bb;
    
    Tmaxv(jjerr) = TT(1);
    
    %% calcolo errore
    errcond(jjerr) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    if jjerr == 1 % soluzione migliore
        TTref = TT;
        drref = dr;
        rrref = rr;
    else
        err(jjerr) = norm(TTref(1:round(dr/drref):end)-TT)/norm(TTref(1:round(dr/drref):end)-Text);
    end
end

%% post production
% temperatura massima
Tmax = TTref(1);
fprintf('Temperatura massima raggiunta dalla pastiglia: %.2f °C.\n',Tmax)

figure
loglog(drv,Tmaxv,'b-.*','linewidth',2)
title('Andamento temperatura massima')
xlabel('\deltar [m]')
ylabel('Errore [-]')
grid on
box on
set(gca,'fontsize',18)

% ditribuzione di temperatura
figure
plot(rrref*1e3,TTref,'b-','linewidth',2)
hold on
plot(rrref(1)*1e3,Tmax,'ko','markersize',10,'markerfacecolor','y')
legend('Andamento temperatura','T_{massima}')
title('Stazionario 1D')
xlabel('Raggio [mm]')
ylabel('Temperatura [°C]')
grid on
box on
set(gca,'fontsize',18)

% andamento errore
figure
loglog(drv,err,'b-.*','linewidth',2)
hold on
loglog(drv,errcond,'r-.*','linewidth',2)
legend('Errore di discretizzazione','Upper bound di troncamento','location','southwest')
title('Errore stazionario 1D')
xlabel('\deltar [m]')
ylabel('Errore [-]')
grid on
box on
set(gca,'fontsize',18)