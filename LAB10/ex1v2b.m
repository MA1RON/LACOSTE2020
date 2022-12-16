clc; clear; close all;
%% dati
% geometrici
DD = 5e-3;
LL = 1; 

% fisici
ro = 140;
cp = 6e3;
Tx0 = 5;
T0t = 10;

uu = .25e-2;

dtv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2];
err = zeros(size(dtv));
errcond = zeros(size(dtv));
for jjerr = 1:length(dtv)
    %% griglia t
    tend = 100;
    dt = dtv(jjerr);
    time = (0:dt:tend)';
    mm = length(time);

    %% griglia x
    dx = LL/100;
    xx = (0:dx:LL)';
    nn = length(xx);

    %% BE - upwind
    TT = ones(nn,1)*Tx0; Tp = TT;

    aa = uu*dt/dx; % numero di Courant CFL
    diag_i = -aa*ones(nn,1);
    diag_p = (1 + aa)*ones(nn,1);
    AA = spdiags([diag_i diag_p],-1:0,nn,nn);

    % --- condizioni al contorno ---
    AA(1,1) = 1;
    AA(1,2) = 0;

    for inst = 2:mm
        bb = Tp;

        % --- condizioni al contorno ---
        bb(1) = T0t;

        TT = AA\bb;
        Tp = TT;
    end
    
    %% valuto errore
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    if jjerr == 1
        dtref = dt;
        TTref = TT;
    else
        err(jjerr) = norm(TTref-TT)/norm(TT-Tx0);
    end
end

%% post production
loglog(dtv,err,'b-.*','linewidth',2)
hold on
loglog(dtv,errcond,'r-.*','linewidth',2)
box on
grid on
set(gca,'fontsize',18)
xlabel('deltat [s]')
ylabel('errore [-]')
title('convergenza temporale pure advection 1D')