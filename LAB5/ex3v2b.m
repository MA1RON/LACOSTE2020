clc; clear; close all;
%% dati
% dominio
% Tint | 0 - ss - contatto - is - end | Text
% geometrici
tss = 1e-3;
tis = 2e-3;

% fisici
Tint = 70;
Text = 20;
hhext = 15;
kkss = 15;
kkis = 1;

dxv = [1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4];
err = zeros(size(dxv));
errcond = zeros(size(dxv));
TTrefnotyetsaved = true;
for jjerr = 1:length(dxv)
    %% griglia
    dx = dxv(jjerr);
    xx = [(0:dx:tss)';(tss+dx:dx:tss+tis)'];
    nncontatto = length(0:dx:tss);
    nn = length(xx);

    %% costruisco il sistema
    diag_si = ones(nn,1);
    diag_p = -2*ones(nn,1);
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
    bb = zeros(nn,1);

    % --- condizioni al contorno ---
    % dirichlet interno
    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = Tint;

    % contatto
    AA(nncontatto,nncontatto-1) = -1;
    AA(nncontatto,nncontatto) = 1+kkis/kkss;
    AA(nncontatto,nncontatto+1) = -kkis/kkss;
    bb(nncontatto) = 0; % pleonastico

    % robin esterno
    AA(end,end-1) = 1;
    AA(end,end) = -1-hhext*dx/kkis;
    bb(end) = -hhext*Text*dx/kkis;

    %% risolvo il sistema
    TT = AA\bb;
    
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    if jjerr == 1 % soluzione più accurata
        TTref = TT;
        xxref = xx;
        dxref = dx;
        nncontattoref = nncontatto;
    else
        err(jjerr) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref(1:round(dx/dxref):end)-Text);
    end
    
    if err(jjerr) < errcond(jjerr) && TTrefnotyetsaved
        TTref = TT;
        xxref = xx;
        dxref = dx;
        nncontattoref = nncontatto;
        TTrefnotyetsaved = false;
    end
end

%% post production
figure('Position', [100 100 1000 700])

subplot(1,2,1)
plot(xxref*1e3,TTref,'b-.','linewidth',2)
hold on
plot([xxref(nncontattoref) xxref(nncontattoref)]*1e3,[min(TTref) max(TTref)],'k--','linewidth',2)
legend('Temperatura', 'Contatto', 'location', 'northeast')
xlabel('Spessore [mm]')
ylabel('Temperatura [°C]')
grid on
box on
set(gca, 'fontsize', 18)
title('Stazionario 1D')

subplot(1,2,2)
loglog(dxv,err,'b-.*','linewidth',2)
hold on
loglog(dxv,errcond,'k--*','linewidth',2)
legend('Discretizzazione', 'Troncamento', 'location', 'southwest')
xlabel('\deltax [m]')
ylabel('Errore [-]')
grid on
box on
set(gca, 'fontsize', 18)
title('Convergenza stazionario 1D')