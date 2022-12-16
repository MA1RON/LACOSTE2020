clc; clear; close all;
%% dati
LL = .1;
alpha = .1e-4;
T0 = 300;
Tx0 = @(x) T0 + 50*sin(pi*x/LL);
Tan = @(x,t) T0 + 50*sin(pi*x/LL).*exp(-pi^2/LL^2*alpha*t);

%% griglia x
dx = .1e-2;
xx = (0:dx:LL)';
nn = length(xx);

%% griglia t
tend = 100;
dt = .1;
time = (0:dt:tend)';
mm = length(time);

ttoshow = 20:20:100;
lineplotted = 0;

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
    
    if time(inst) >= ttoshow(lineplotted + 1)
        plot(xx,TT,'-','Color',[1-inst/mm 0 inst/mm],'linewidth',2)
        hold on
        
        lineplotted = lineplotted +1;
    end
    
%     plot(xx,TT,'-','Color',[1-inst/mm 0 inst/mm],'linewidth',2)
%     xlim([0 .1])
%     ylim([300 400])
%     drawnow
    
    Tp = TT;
end