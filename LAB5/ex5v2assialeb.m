clc; % clear; close all;
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

%% ciclo per errori
dxv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2];
err = zeros(size(dxv));
errcond = zeros(size(dxv));
for jjerr = 1:length(dxv)
    %% griglia spaziale
    dx = dxv(jjerr);
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
    
    %% controllo errori
    errcond(jjerr) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    if jjerr == 1 % soluzione più attendibile
        TTref = TT;
        dxref = dx;
        xxref = xx;
    else
        err(jjerr) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref-Tinf);
    end
end

%% post production
figure('Position', [10 10 1300 600])
subplot(1,2,1)
plot(xxref*1e2,TTref,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Lunghezza (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)

subplot(1,2,2)
loglog(dxv,err,'b-.*','linewidth',2)
hold on
loglog(dxv,errcond,'r-.*','linewidth',2)
legend('Errore discretizzazione','Upper bound troncamento','location','northeastoutside')
grid on
box on
title ('Andamento errore')
xlabel('\Deltax (m)')
ylabel('Errore (-)')
set(gca,'fontsize',18)