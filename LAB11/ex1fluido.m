clc; clear; close all;
%% dati
% geometrici
DDint = 1e-2;
DDout = DDint+2*1e-3;
LL = 2;

surfint = pi*DDint^2/4;

% fisici
GG = .12;
T0 = 25;
hh = 2e3;

roac = 7800;
cpac = 500;
kkac = 30;

rowa = 997;
cpwa = 4186;

uu = GG/surfint/rowa;

q0 = 2e3;
qqq = @(x) q0*exp(-5*(x-LL/2).^2)/surfint;

%% griglia temporale
timetoshow = [1 5 10]; lineplotted = 0;
tend = timetoshow(end);
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% BE - upwind
TT = ones(nn,1)*T0; Tp = TT;

aa = uu*dt/dx;
diag_i = -aa*ones(nn,1);
diag_p = (1+aa)*ones(nn,1);
AA = spdiags([diag_i diag_p],-1:0,nn,nn);

% --- condizioni al contorno ---
AA(1,1) = 1;
AA(1,2) = 0;

for inst = 2:mm
    bb = Tp + qqq(xx)/rowa/cpwa*dt;
    
    % --- condizioni al contorno ---
    bb(1) = T0;
    
    % --- risolvo il sistema ---
    TT = AA\bb;
    Tp = TT;
    
    figure(1)
    plot(xx,TT,'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
    title('Transitorio 1D')
    xlabel('x (m)')
    ylabel('Temperature (^oC)')
    grid on
    box on
    set(gca,'fontsize',18)
    ylim([T0 28.5])
    xlim([0 LL])
    drawnow
    
    %{
    if time(inst) >= timetoshow(lineplotted + 1)
        figure(2)
        plot(xx,TT,'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
        hold on
        title('Transitorio 1D')
        xlabel('x (m)')
        ylabel('Temperature (^oC)')
        grid on
        box on
        set(gca,'fontsize',18)
        lineplotted = lineplotted +1;
    end
    %}
end