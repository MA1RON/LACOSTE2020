clc; clear; close all;

%% dati
delta = 50e-3; %m
kk = 5; %W/mK
qqq = @(x) 5e5*sin(pi*x/delta); %W/m^3
Ta = 300; %K
Tb = 300; %K

%% analitica
C2 = Ta;
C1 = (Tb-C2)/delta;
Tan = @(x) 5e5/kk*sin(pi/delta*x)*delta^2/pi^2 + C1*x + C2;

%% ciclo errori
dxv = [1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3]; %in ordine crescente!

%% errori
toll = .1e-2;
err_an = zeros(size(dxv));
err_num = zeros(size(dxv));
RR_an = zeros(size(dxv));
RR_num = zeros(size(dxv));
for jj = 1:length(dxv)
    %% griglia
    dx = dxv(jj);
    xx = (0:dx:delta)';
    nn = length(xx);
    
    %% AA
    my_diag = ones(nn,1);
    AA = spdiags([my_diag -2*my_diag my_diag], -1:1, nn, nn);
    bb = qqq(xx)/-kk*dx^2;

    %% coco
    AA(1,1) = 1; AA(1,2) = 0; bb(1) = Ta; %dirichlet
    AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = Tb; %dirichlet

    %% TT & err
    TT = AA\bb;
    
    err_an(jj) = norm((TT-Ta)-(Tan(xx)-Ta))/norm(Tan(xx)-Ta);
    RR_an(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    
    if jj == 1
        dxref = dx;
        Tref = TT;
        
        err_num(jj) = 0;
    else
        err_num(jj) = norm(TT-Tref(1:round(dx/dxref):end))/norm(Tref(1:round(dx/dxref):end)-Ta);
    end
end

%% plotto soluzione
%{
plot(xx, TT, 'b-', 'linewidth', 2)
hold on 
plot(xx, Tan(xx), 'r-.', 'linewidth', 2)
legend('numerica', 'analitica')
grid on
xlabel('x [m]')
ylabel('T [K]')
set(gca, 'fontsize', 18)
%}

%% plotto errore

loglog(dxv, err_an, '-o', 'linewidth', 2)
hold on 
loglog(dxv, err_num, '-s', 'linewidth', 2)
hold on
loglog(dxv, RR_an, 'r-.', 'linewidth', 2)
legend('discretizzazione analitica', 'discretizzazione numerica', 'limite sup arrotondamento', 'location', 'south')
grid on
xlabel('x [m]')
ylabel('T [K]')
set(gca, 'fontsize', 18)
%}