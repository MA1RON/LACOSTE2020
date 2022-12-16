clc; clear; close all;
%% dati
% geometrici
DD = 5e-3;
LL = 1;

% fisici
Tx0 = 5;
T0t = 10;

uu = @(T) 10./(T+1);

%% griglia t
timetoshow = .25:.25:1;
lineplotted = 0;
tend = timetoshow(end);
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% BE - upwind
TT = ones(nn,1)*Tx0; Tp = TT;

for inst = 2:mm
    % frozen coefficients
    aa = uu(Tp)*dt/dx; % numero di Courant CFL
    diag_i = [-aa(2:end);0];
    diag_p = (1 + aa);
    AA = spdiags([diag_i diag_p],-1:0,nn,nn);
    bb = Tp;

    % --- condizioni al contorno ---
    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = T0t;
    
    TT = AA\bb;
    Tp = TT;
    
    if time(inst) >= timetoshow(lineplotted +1)
        plot(xx,TT,'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
        box on
        grid on
        set(gca,'fontsize',18)
        xlabel('x [m]')
        ylabel('T [K]')
        title('Pure advection 1D')
        
        hold on
        
        lineplotted = lineplotted +1;
    end
    %{
    plot(xx,TT,'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
    box on
    grid on
    set(gca,'fontsize',18)
    xlabel('x [m]')
    ylabel('T [K]')
    title('Pure advection 1D')
    ylim([-5 20])
    xlim([0 LL])
    drawnow
    %}
end