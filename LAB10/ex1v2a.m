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

%% griglia t
timetoshow = 25:25:100;
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
    plot(xx,TT)
    drawnow
    %}
end
