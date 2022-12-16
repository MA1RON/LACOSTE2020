clc; clear; close all;

%% dati
kk = .55; %W/mK
tb = .25; %m
Ti = 20; %°C
Te = -5; %°C
hh = 25; %W/m^2K
qqq = 1e3; %W/m^3

% to find: Tmax, x_Tmax, T_xmedio, errori.

%% analitica
C2 = Ti;
C1 = (qqq*tb/kk+hh*qqq*tb^2/2/kk^2-C2*hh/kk+hh*Te/kk)/(1+hh*tb/kk);
Tan = @(x) qqq/-kk*x.^2/2+C1*x+C2;

%% griglie
dxv = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3];
err_num_an = ones(size(dxv));
err_num_an_max = ones(size(dxv));
err_num_an_med = ones(size(dxv));
err_num_ref = ones(size(dxv));
RR = ones(size(dxv));
for jj = 1:length(dxv)
    dx = dxv(jj);
    xx = (0:dx:tb)';
    nn = length(xx);

    %% AA & bb
    my_diag = ones(nn,1);
    AA = spdiags([my_diag -2*my_diag my_diag], -1:1, nn, nn);
    bb = qqq/-kk*dx^2*ones(nn,1);

    %% coco & TT
    AA(1,1) = 1; AA(1,2) = 0; bb(1) = Ti; %dirichlet
    AA(end,end) = 1+hh*dx/kk; AA(end,end-1) = -1; bb(end) = hh*dx*Te/kk; %robin

    TT = AA\bb;
    Tmax = max(TT);
    Tmed = TT(round(nn/2));
    
    %% errori
    err_num_an(jj) = norm(TT-Tan(xx))/norm(Tan(xx)-Ti);
    err_num_an_max(jj) = norm(Tmax-max(Tan(xx)))/norm(max(Tan(xx))-Ti);
    err_num_an_med(jj) = norm(Tmed-Tan(round(xx/2)))/norm(Tan(round(xx/2))-Ti);
    
    if jj == 1
        Tref = TT;
        dxref = dx;
    else
        err_num_ref(jj) = norm(TT-Tref(1:round(dx/dxref):end))/norm(Tref(1:round(dx/dxref):end)-Ti);
    end
    
    RR(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
end

%% plotto soluzione
plot(xx,TT,'b-','linewidth',2)
hold on
plot(xx,Tan(xx),'r-.','linewidth',2)
legend('soluzione numerica', 'soluzione analitica')
grid on
xlabel('x [m]')
ylabel('T [°C]')

%% plotto errori
figure()
loglog(dxv,err_num_an,'b-s','linewidth',2)
hold on
loglog(dxv(2:end), err_num_ref(2:end), 'k-s','linewidth',2)
hold on
loglog(dxv, RR, 'r-.s','linewidth',2)
hold on
loglog(dxv,err_num_an_max,'g-s','linewidth',2)
hold on
loglog(dxv,err_num_an_med,'c-s','linewidth',2)
legend('errore numerico analitica', ...
    'errore numerico riferimento',...
    'limite di condizionamento (arrotondamento)',...
    'errore sul massimo',...
    'errore sul medio',...
    'location','south')
grid on
xlabel('dx [m]')
ylabel('err [1]')