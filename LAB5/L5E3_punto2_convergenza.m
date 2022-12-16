clc; clear; close all;

%% dati
tu = 2e-3; %m
td = 1e-3; %m
Tu = 20; %°C
Td = 70; %°C
hh = 15; %W/m^2K
ku = 1; %W/mK
kd = 15; %W/mK

%% ciclo errori
dxv=[1.5e-7 1.5e-6 1.5e-5 1.5e-4];

err_ref = ones(length(dxv),1);
RR = ones(length(dxv),1);
%% griglia {(up) 0->robin->i->s->dirichlet->end (down)}
for jj = 1:length(dxv)
    dx = dxv(jj); %m
    xu = (0:dx:tu)';
    nu = length(xu);
    xd = (tu+dx:dx:tu+td)';
    xx = [xu;xd];
    nn = length(xx);

    %% AA & bb
    my_diag = ones(nn,1);
    AA = spdiags([my_diag -2*my_diag my_diag], -1:1, nn, nn);
    bb = zeros(nn,1);

    %% coco & TT
    %robin
    AA(1,1) = hh/ku*dx+1;
    AA(1,2) = -1;
    bb(1) = hh*Tu*dx/ku;

    %continuità
    AA(nu,nu-1) = -ku;
    AA(nu,nu) = ku+kd;
    AA(nu,nu+1) = -kd;
    bb(nu) = 0;

    %dirichlet
    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = Td;

    TT = AA\bb;
    
    %% errori
    RR (jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    if jj == 1
        TTref = TT;
        dxref = dx;
        xxref = xx;
    else
        err_ref(jj) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref(1:round(dx/dxref):end)-Td);
    end
end

%% plotto soluzione
plot(xxref,TTref, 'b-', 'linewidth', 2)
hold on
plot([2e-3 2e-3], [70 68.5], 'k-.', 'linewidth', 2)
legend('Temperatura', 'Contatto', 'location', 'south')
xlabel('x [m]')
ylabel('T [°C]')
grid on
axis([0 3e-3 68.5 70])
set(gca, 'fontsize', 18)

close

%% plotto errore
figure()
loglog(dxv(2:end),err_ref(2:end), 'b-s', 'linewidth', 2)
hold on
loglog(dxv,RR, 'r-.s', 'linewidth', 2)
legend('troncamento', 'arrontamento')
xlabel('dx [m]')
ylabel('err [1]')
grid on