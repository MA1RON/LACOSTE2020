clc; clear; close all;

%% dati
ri = 5e-3; %m
re = 1e-2; %m
kk = .1; %W/mK
qqq = 2e6; %W/m^3
Tw = 300; %°C
hh = 10; %W/m^2K

%% analitica
C1 = qqq*ri^2/2/kk;
T2 = (qqq*re/-2/kk+C1/re)/hh*-kk+Tw;
C2 = T2+qqq*re^2/kk/4-C1*log(re);

Tan = @(r) qqq*r.^2/-kk/4+C1*log(r)+C2;

%% ciclo errori
drv = [5e-7 1e-6 2e-6 5e-6 1e-5 2e-5];
err_num_an = ones(size(drv));
err_num_ref = ones(size(drv));
err_num_max = ones(size(drv));
RR = ones(size(drv));
for jj = 1:length(drv)
    %% griglia
    dr = drv(jj);
    rr = (ri:dr:re)';
    nn = length(rr);

    %% AA & bb
    diag_s = dr+2*rr;
    diag_p = -4*rr;
    diag_i = -dr+2*rr;
    AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);
    bb = -2*qqq/kk*dr^2*rr;

    %% coco % TT
    AA(1,1) = 1; AA(1,2) = -1; bb(1) = 0; %neumann omo
    AA(end,end-1) = 1; AA(end,end) = hh*dr/-kk-1; bb(end) = Tw*hh*dr/-kk; %robin
    TT = AA\bb;
    Tmax = max(TT);
    
    %% errori
    err_num_an(jj) = norm(TT-Tan(rr))/norm(Tan(rr)-Tw);
    err_num_max(jj) = norm(Tmax-max(Tan(rr)))/norm(max(Tan(rr))-Tw);
    if jj == 1
        Tref = TT;
        drref = dr;
    else
        err_num_ref(jj) = norm(TT-Tref(1:round(dr/drref):end))/norm(Tref(1:round(dr/drref))-Tw);
    end
    RR(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
end

%% plotto soluzione
plot(rr,Tan(rr),'r-.','linewidth', 2)
hold on
plot(rr, TT,'b-','linewidth', 2)
legend('soluzione analitica', 'soluzione numerica')
grid on
xlabel('x [m]')
ylabel('T [°C]')
set(gca, 'fontsize', 18)

%% plotto errori
figure()
loglog(drv,err_num_an,'r-s','linewidth', 2)
hold on
loglog(drv,err_num_max,'b-s','linewidth', 2)
hold on
loglog(drv(2:end),err_num_ref(2:end),'k-s','linewidth', 2)
hold on
loglog(drv,RR,'g-s','linewidth', 2)
legend('errore su analitica',...
    'errore del massimo',...
    'errore su riferimento',...
    'errore sul condizionamento',...
    'location','south')
grid on
xlabel('dx [m]')
ylabel('err [1]')
set(gca, 'fontsize', 18)