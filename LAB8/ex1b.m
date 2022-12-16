clc; clear; close all;
%% dati
LL = .1;
alpha = .1e-4;
T0 = 300;
Tx0 = @(x) T0 + 50*sin(pi*x/LL);
Tan = @(x,t) T0 + 50*sin(pi*x/LL).*exp(-pi^2/LL^2*alpha*t);

%% griglia x
dx = .01;
xx = (0:dx:LL)';
nn = length(xx);

%% ciclo per errori
dtv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1];
err = zeros(size(dtv));
errcond = zeros(size(dtv));
for jjerr = 1:length(dtv)
    %% griglia t
    tend = 100;
    dt = dtv(jjerr);
    time = (0:dt:tend)';
    mm = length(time);

    %% BE
    aa = dt*alpha/dx^2;
    diag_si = -aa*ones(nn,1);
    diag_p = (1+2*aa)*ones(nn,1);
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);

    % --- condizioni al contorno ---
    % dirichlet
    AA(1,1) = 1;
    AA(1,2) = 0;

    % dirichlet
    AA(end,end-1) = 0;
    AA(end,end) = 1;

    TT = Tx0(xx); Tp = TT;
    for inst = 2:mm
        bb = Tp;

        % --- condizioni al contorno ---
        % dirichlet
        bb(1) = T0;

        % dirichlet
        bb(end) = T0;

        TT = AA\bb;

        Tp = TT;
    end
    
    %% valuto errore
    errcond(jjerr) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    if jjerr == 1
        TTref = TT;
    else
        err(jjerr) = norm(TTref-TT)/norm(TTref-T0);
    end
end

loglog(dtv,errcond,'b-.*','linewidth',2)
hold on
loglog(dtv,err,'r-.*','linewidth',2)