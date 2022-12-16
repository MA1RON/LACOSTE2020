clc; clear; close all;
%% dati
% geometrici
DDint = 2e-3;
DDout = 4e-3;

surf = pi*(DDout^2-DDint^2)/4;

% fisici
II = 10; JJ = II/surf;
roel = 2.5e-7;
alpha = 11e-5;
kk = 3.5;
rocp = kk/alpha;
T0 = 300;
hh = 100;

qqq = roel*JJ^2;
tol = 1e-5;

%% griglia t
tend = 17; % guess
dt = tend/50;
time = (0:dt:tend)';
mm = length(time);

%% griglia x
dr = (DDout-DDint)/200;
rr = (DDint/2:dr:DDout)';
nn = length(rr);

%% BE
timetoplot = linspace(time(1),time(end),10);
lineplotted = 0;

TT = ones(nn,1)*T0;
Tp = TT;
aa = dt*alpha;
diag_s = -aa/dr^2-aa/2./rr/dr;
diag_p = (1+2*aa/dr^2)*ones(nn,1);
diag_i = -aa/dr^2+aa/2./rr/dr;
AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1,nn,nn);

% --- condizioni al contorno ---
% robin
AA(1,1) = -1-hh/kk*dr;
AA(1,2) = 1;

% neumann
AA(end,end-1) = 1;
AA(end,end) = -1;

inst = 2; res = 1;
while inst<mm && res > tol
    bb = Tp + qqq/rocp*dt;
    
    % --- condizioni al contorno ---
    % robin
    bb(1) = hh/kk*dr*-T0;
    % neumann
    bb(end) = 0;
    
    TT = AA\bb;
    
    %{
    plot(rr,TT)
    ylim([300 330])
    drawnow
    %}
    %{
    if time(inst) >= timetoplot(lineplotted + 1)
        plot(rr*1e3,TT,'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
        xlabel('Raggio [mm]')
        ylabel('Temperatura [°C]')
        grid on
        box on
        set(gca,'fontsize',18)
        ylim([300 330])
        
        hold on
        
        lineplotted = lineplotted +1;
    end
    %}
    res = norm(TT-Tp)/norm(TT-T0);
    Tp = TT;
    inst = inst + 1;
end
