clc; clear; close all;
%% dati
LL = 50e-3;
kk = 5;
q0 = 500e3;
qqq = @(x) q0*sin(pi*x/LL);
Tside = 300;

tol = .1e-2;
tol_raggiunta = false;

%% soluzione analitica
Tan = @(x) q0/kk*LL^2/pi^2*sin(pi*x/LL) + Tside;

%% controllo la tolleranza
dxv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3];
err = zeros(size(dxv));
errcond = zeros(size(dxv));
errnum = zeros(size(dxv));
for jjerr = 1:length(dxv)
    %% griglia spaziale
    dx = dxv(jjerr);
    xx = (0:dx:LL)';
    nn = length(xx);

    %% creo il problema
    diag_si = ones(nn,1);
    diag_p = -2*ones(nn,1);
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
    bb = - qqq(xx)/kk*dx^2;

    %% condizioni al contorno
    AA(1,1) = 1; AA(1,2) = 0; bb(1) = Tside;
    AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = Tside;

    %% soluzione numerica
    TT = AA\bb;
    
    %% controllo l'errore
    err(jjerr) = norm(TT-Tan(xx))/norm(Tan(xx)-Tside);
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    
    if jjerr == 1 % soluzione più accurata
        xxnumref = xx;
        TTnumref = TT;
        dxnumref = dx;
    else
        errnum(jjerr) = norm(TTnumref(1:round(dx/dxnumref):end)-TT)/norm(TTnumref(1:round(dx/dxnumref):end)-Tside);
    end
    
    if err(jjerr) <= tol && not(tol_raggiunta)
        xxref = xx;
        TTref = TT;
        errref = err(jjerr);
        tol_raggiunta = true;
    end
end
if tol_raggiunta
    fprintf('Tolleranza raggiunta:\nErrore = %.2e.\nTolleranza = %.2e.\n',errref,tol)
else
    fprintf('Tolleranza non raggiunta.\n')
    xxref = xx;
    TTref = TT;
end

%% post production
% grafico temperatura
plot(xxref*1e3,TTref-273.15,'b-','linewidth',2)
hold on
plot(xxref*1e3,Tan(xxref)-273.15,'r-.','linewidth',2)
legend('Soluzione numerica','Soluzione analitica','location','south')
grid on
box on
xlabel('Spazio [mm]')
ylabel('Temperatura [°C]')
set(gca,'fontsize',18)

% grafico errore
figure
loglog(dxv,err,'b*','linewidth',2)
hold on
loglog(dxv,errnum,'g*','linewidth',2)
hold on
loglog(dxv,errcond,'r*','linewidth',2)
legend('Errore sulla analitica','Errore sulla numerica','Errore di troncamento','location','northwest')
grid on
box on
xlabel('\deltax [m]')
ylabel('Errore [-]')
set(gca,'fontsize',18)

%% conservazione dell'energia
Ein = trapz(xx,qqq(xx))
Eout = 2*kk*(TT(end-1)-TT(end))/dx