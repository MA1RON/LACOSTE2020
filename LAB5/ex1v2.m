clc; clear; close all;
%% dati
% dominio
% Text | 0 - is - contatto - br - end | Tint
% geometrici
tb = .25;
ti = .1;

% fisici
kb = .55;
Tint = 20;
Text = -5;
hhext = 25;
kiref = .04;

%% griglia spaziale
dx = 1e-3;
xx = [(0:dx:ti)';(ti+dx:dx:tb)'];
ncont = length(0:dx:ti);
nn = length(xx);

%% creo il sistema
diag_si = ones(nn,1);
diag_p = -2*ones(nn,1);
AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
bb = zeros(nn,1);

%% ciclo sulla conduttività dell'isolante
kiv = sort([linspace(.004,.55), kiref]);
qqlost = zeros(size(kiv));
for jjki = 1:length(kiv)
    ki = kiv(jjki);

    % condizioni al contorno
    % robin esterno
    AA(1,1) = 1+dx*hhext/ki;
    AA(1,2) = -1;
    bb(1) = dx*hhext/ki*Text;

    % contatto
    AA(ncont,ncont-1) = -1;
    AA(ncont,ncont) = 1+kb/ki;
    AA(ncont,ncont+1) = -kb/ki;
    % bb(cont) = 0;

    % dirichlet interno
    AA(end,end-1) = 0;
    AA(end) = 1;
    bb(end) = Tint;

    %% risolvo il sistema
    TT = AA\bb;
    
    %% calcolo il calore disperso
    qqlost(jjki) = abs(hhext*(Text-TT(1))); % W/m^2
    if ki == kiref
        TTref = TT;
        qqlostref = qqlost(jjki);
    end
end

%% posto production
% andamento temperatura
figure
plot(xx*1e2,TTref,'b-.','linewidth',2)
hold on
plot([ti ti]*1e2,[min(TT) max(TT)],'k--','linewidth',2)
legend('Temperatura','Superficie di contatto','location','southeast')
grid on
box on
xlabel('Spessore [cm]')
ylabel('Temperatura [^oC]')
title('Stazionario 1D')
set(gca,'fontsize',16)

% andamento flusso disperso
figure
plot(kiv,qqlost,'b-.','linewidth',2)
hold on
plot(kiref,qqlostref,'ko','markersize',10,'markerfacecolor','y')
grid on
box on
xlabel('Conducibilità isolante [W/m^2K]')
ylabel('Flusso disperso [W/m^2]')
title('Flusso perso stazionario 1D')
set(gca,'fontsize',16)
