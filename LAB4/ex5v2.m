clc; clear; close all;
%% dati
% geometrici
DDint = 20e-2;
DDext = 30e-2;

% fisici
Text = 5;
hhext = 100;
qqint = 50;
kk = .2;
qqq = 0;

%% soluzione analitica
C1 = -qqint*(DDint/2)^2/kk;
T2 = -kk/hhext*C1/(DDext/2)^2+Text;
C2 = T2+C1/(DDext/2);
Tan = @(r) -C1./r + C2;

%% ciclo per l'errore
drv = [5e-5 1e-4 2e-4 5e-4 1e-3];
erran = zeros(size(drv));
errnum = zeros(size(drv));
errcond = zeros(size(drv));
flux = zeros(size(drv));
for jjerr = 1:length(drv)
    %% griglia
    dr = drv(jjerr);
    rr = (DDint/2:dr:DDext/2)';
    nn = length(rr);

    %% costruisco il sistema
    diag_s = 1+rr/dr;
    diag_p = -2*rr/dr;
    diag_i = -1+rr/dr;
    AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1,nn,nn);
    bb = rr*dr*-qqq/kk;

    % --- condizioni al contorno ---
    % neumann interno
    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = -qqint*dr/kk;

    % robin esterno
    AA(end,end-1) = 1;
    AA(end,end) = -1-hhext*dr/kk;
    bb(end) = -hhext*dr*Text/kk;

    %% risolvo il sistema
    TT = AA\bb;
    
    flux(jjerr) = hhext*(TT(end)-Text);
    
    %% calcolo errore
    errcond(jjerr) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    erran(jjerr) = norm(TT-Tan(rr))/norm(Tan(rr)-Text);
    if jjerr == 1 % soluzione più precisa
        drref = dr;
        TTref = TT;
        rrref = rr;
        fprintf('Flusso scambiato sulla superficie esterna: %.2f W/m^2.\n',flux(jjerr))
        fprintf('Potenza totale dispersa: %.2f W.\n',flux(jjerr)*pi*DDext^2)
    else
        errnum(jjerr) = norm(TT-TTref(1:round(dr/drref):end))/norm(TTref(1:round(dr/drref):end)-Text);
    end
end

%% post production
% andamento flusso
figure
loglog(drv, flux,'r-.*','linewidth', 2)
grid on
box on
title('Andamento flusso di calore')
xlabel('\deltar [m]')
ylabel('Flusso [W/m^2]')
set(gca, 'fontsize', 18)

% andamento temperatura
figure
plot(rrref*1e2, TTref,'b-','linewidth', 2)
hold on
plot(rrref*1e2,Tan(rrref),'r-.','linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
box on
title('Stazionario 1D')
xlabel('Raggio [cm]')
ylabel('Temperatura [°C]')
set(gca, 'fontsize', 18)

% andamento errore
figure
loglog(drv, errcond,'r-.*','linewidth', 2)
hold on
loglog(drv, erran,'b-.*','linewidth', 2)
hold on
loglog(drv, errnum,'g-.*','linewidth', 2)
legend('Upper bound di condizionamento', "Errore sull'analitica",'Errore sulla numerica',...
    'location','northeastoutside')
grid on
box on
title('Errore stazionario 1D')
xlabel('\deltar [m]')
ylabel('Errore [-]')
set(gca, 'fontsize', 18)
