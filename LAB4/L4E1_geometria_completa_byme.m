clc; clear; close all;

%% dati
delta = 50e-3; %m
kk = 5; %W/mK
qqq = @(x) 5e5*sin(pi*x/delta); %W/m^3
Ta = 300; %K
Tb = 300; %K

toll = .1e-2;
n_err_max = 7;
n_err_min = 3;
err = ones(3*(n_err_max-n_err_min), 1);
RR = ones(3*(n_err_max-n_err_min), 1);

%% analitica
C2 = Ta;
C1 = (Tb-C2)/delta;
Tan = @(x) 5e5/kk*sin(pi/delta*x)*delta^2/pi^2 + C1*x + C2;

%% griglie
jj = 1;
dxv = logspace(-n_err_max,-n_err_min,n_err_max-n_err_min);
dxv = sort([dxv 2*dxv 5*dxv],'descend');
for dx = dxv
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
    
    err(jj) = norm(TT-Tan(xx))/norm(Tan(xx)-Ta);
    RR(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    jj = jj+1;
end

%% plotto soluzione

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
%{
loglog(dxv, err, 'b-', 'linewidth', 2)
hold on 
plot(dxv, RR, 'r-.', 'linewidth', 2)
legend('discretizzazione', 'limite sup arrotondamento', 'location', 'south')
grid on
xlabel('x [m]')
ylabel('T [K]')
set(gca, 'fontsize', 18)
%}