clc; clear; close all;
%% dati
% dominio
% Text:0 --- muro --- end:Tint
% geometrici
tb = 25e-2;

% fisici
kb = .55;
Tint = 20;
hhext = 25;
Text = -5;
qqq = 1e3;

%% verifico convergenza dell'errore
dxv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2];
err = zeros(size(dxv));
err2 = zeros(size(dxv));
errcond = zeros(size(dxv));
for jjerr = 1:length(dxv)
    %% griglia
    dx = dxv(jjerr);
    xx = (0:dx:tb)';
    nn = length(xx);

    %% costruisco il sistema
    diag_si = ones(nn,1);
    diag_p = -2*ones(nn,1);
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
    bb = -qqq/kb*ones(nn,1)*dx^2;

    % --- condizioni al contorno ---
    % robin esterno
    AA(1,1) = -1+hhext*dx;
    AA(1,2) = 1;
    bb(1) = hhext*dx*Text;
    
    % dirichlet interno
    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = Tint;
    
    % robin esterno polinomiale
    AA2 = AA; bb2 = bb;
    AA2(1,1) = -1+hhext*dx;
    AA2(1,2) = 4;
    AA2(1,3) = -3;
    bb2(1) = 2*hhext*dx*Text;
    
    %% risolvo il sistema
    TT = AA\bb;
    TT2 = AA2\bb2;
    
    %% calcolo errori
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    if jjerr == 1 % soluzione più precisa
        dxref = dx;
        TTref = TT;
        TT2ref = TT2;
        xxref = xx;
    else
        err(jjerr) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref(1:round(dx/dxref):end)-Tint);
        err2(jjerr) = norm(TT2ref(1:round(dx/dxref):end)-TT2)/norm(TT2ref(1:round(dx/dxref):end)-Tint);
    end
end

%% post production
% temperatura massima
[Tmax, nnmax] = max(TTref);
xmax = tb-xxref(nnmax);
fprintf('La temperatura massima raggiunta è %.2f °C, a %.2f cm dalla parete interna.\n',Tmax,xmax*1e2)

% temperatura centrale
xxave = xxref(round(length(xxref)/2));
TTave = TTref(round(length(xxref)/2));
fprintf('La temperatura del punto medio è %.2f °C.\n',TTave)

% plot temperatura
figure
plot(xxref*1e2,TTref,'b-.','linewidth',2)
hold on
plot(xxref(nnmax)*1e2,TTref(nnmax),'ko','markersize',10,'markerfacecolor','y')
hold on
plot(xxave*1e2,TTave,'ko','markersize',10,'markerfacecolor','g')
legend('Distribuzione di temperatura','T_{max}','T_{punto medio}','location','southeast')
xlabel('Lunghezza [cm]')
ylabel('Temperatura [°C]')
title('Stazionario 1D')
grid on
box on
set(gca,'fontsize',18)

% plot errore
figure
loglog(dxv,err,'b-.*','linewidth',2)
hold on
loglog(dxv,err2,'g-.*','linewidth',2)
hold on
loglog(dxv,errcond,'r-.*','linewidth',2)
legend('In avanti','Polinomiale','Upper bound di condizionamento','location','northeastoutside')
xlabel('\deltax [m]')
ylabel('Errore [-]')
title('Errore stazionario 1D')
grid on
box on
set(gca,'fontsize',18)