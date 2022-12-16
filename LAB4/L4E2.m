clc; clear; close all;

%% dati
DD = .01; %m
LL = 4; %m
II = 1e3; %A
Tb = 77; %K
hh = 500; %W/m^2K
roel = 1.75e-8; %hom*m
kk = 350; %W/mK

SE = pi*DD*LL;
SI = pi*DD^2/4;
VV = SI*LL;

RR = roel*LL/(pi*DD^2/4); %hom
qqq = (RR*II^2)/VV + SE*hh/VV*Tb; %W/m^3

%% griglia
dxv = [1e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2];
err_an = ones(length(dxv),1);
err_num = ones(length(dxv),1);
RA = ones(length(dxv),1);
for jj = 1:length(dxv)
    dx = dxv(jj); 
    xx = LL/2:dx:LL;
    nn = length(xx);

    %% AA & bb
    s_diag = ones(nn,1);
    c_diag = (-2-SE*hh/kk/VV*dx^2)*ones(nn,1);
    i_diag = ones(nn,1);

    AA = spdiags([s_diag c_diag i_diag], -1:1, nn, nn);
    bb = -qqq*dx^2/kk*ones(nn,1);

    %% coco
    AA(1,1) = -1; AA(1,2) = 1; bb(1) = 0; % neumann omo simm
    AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = Tb;

    TT = AA\bb;
    
    %% err
    if jj == 1
        dxref = dx;
        Tref = TT;
        
        err_num(jj) = 0;
    else
        RA(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
        err_num(jj) = norm(TT-Tref(1:round(dx/dxref):end))/norm(Tref(1:round(dx/dxref):end)-Tb);
    end
end

%% plotto soluzione

plot(xx,TT,'b-','linewidth', 2)
xlabel('x [m]')
ylabel('T [K]')
grid on
set(gca, 'fontsize', 18)
%}

%% plotto errore
figure()
loglog(dxv,err_num,'b-s','linewidth', 2)
hold on
loglog(dxv, RA,'r-.s','linewidth', 2)
xlabel('dx [m]')
ylabel('err [1]')
grid on
%axis([LL/2, LL, 76, 92])
%}