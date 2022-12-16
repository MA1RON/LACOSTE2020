clc; clear; close all;

%% dati
Ta = 5; %°C
hh = 100; %W/m^2K
qq = 50; %W/m^2
ri = .1; %m
re = .15; %m
kk = .2; %W/mK

%% analitica
C1 = -qq*ri^2/kk;
T2 = -kk/hh*C1/re^2+Ta;
C2 = T2+C1/re;
Tan = @(r) -C1./r + C2;

%% griglia
drv = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3];
err_num_an = ones(size(drv));
err_num_ref = ones(size(drv));
RR = ones(size(drv));
for jj = 1:length(drv)
    dr = drv(jj); 
    rr = (ri:dr:re)';
    nn = length(rr);

    %% AA & bb
    diag_i = dr+rr;
    diag_p = -2*rr;
    diag_s = -dr+rr;
    AA = spdiags([[diag_s(2:end); 0] diag_p [0;diag_i(1:end-1)]], -1:1, nn, nn);
    % errato prima diag_i poi diag_s
    bb = zeros(nn,1);

    %% coco & TT
    %neumann non omo
    AA(1,1) = -1; 
    AA(1,2) = 1; 
    bb(1) = qq/-kk*dr; 
    %robin
    AA(end,end-1) = -1; 
    AA(end,end) = 1 + hh/kk*dr; 
    bb(end) = Ta*hh*dr/kk; 

    TT = AA\bb;
    
    %% errori
    err_num_an(jj) = norm(TT-Tan(rr))/norm(Tan(rr)-Ta);
    if jj == 1
        Tref = TT;
        drref = dr;
    else
        err_num_ref(jj) = norm(TT-Tref(1:round(dr/drref):end))/norm(Tref(1:round(dr/drref):end)-Ta);
    end
    RR(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
end

%% flusso in re (richiesta)
FE = -kk*(TT(end)-TT(end-1))/dr; % W/m^2

%% plotto soluzione
plot(rr, TT,'b-','linewidth', 2)
hold on
plot(rr,Tan(rr),'r-.','linewidth', 2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('x [m]')
ylabel('T [°C]')
set(gca, 'fontsize', 18)

%% plotto errori
figure()
loglog(drv, err_num_an,'b-s','linewidth', 2)
hold on
loglog(drv(2:end),err_num_ref(2:end),'r-s','linewidth', 2)
hold on
loglog(drv,RR,'k-.s','linewidth', 2)
legend('errore numerico analitica', 'errore numerico ref','limite di condizionamento',...
    'location', 'south')
grid on
xlabel('dx [m]')
ylabel('err [1]')
set(gca, 'fontsize', 18)