clc; clear; close all;

%% dati
tb = .25; %m
ti = .1; %m
kb = .55; %W/mK
ki = .04; %W/mK
Ti = 20; %°C
Te = -5; %°C
hh = 25; %W/m^2K

%% griglie
dxi = 1e-3;
dxb = 1e-3;
xi = (0:dxi:ti)';
xb = (ti+dxb:dxb:ti+tb)';
xx = [xi;xb];
ni = length(xi);
nb = length(xb);
nn = length(xx); % ni+nb

%% perdite
dk = 1e-2;
dkv = (ki:dk:kb);
nnk = length(dkv);

qqp = ones(nnk,1);

for jj = 1:nnk
    %% AA & bb
    diag_s = ones(nn,1);
    diag_p = -2*ones(nn,1);
    diag_i = ones(nn,1);

    AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);

    bb = zeros(nn,1);

    %% coco & TT
    % robin
    AA(1,1) = hh/dkv(jj)*dxi+1;
    AA(1,2) = -1;
    bb(1) = Te*hh/dkv(jj)*dxi;

    % flusso continuo
    AA(ni,ni-1) = -dkv(jj)/dxi;
    AA(ni,ni) = dkv(jj)/dxi+kb/dxb;
    AA(ni, ni+1) = -kb/dxb;
    bb(ni) = 0;

    % dirichlet
    AA(end, end-1) = 0;
    AA(end, end) = 1;
    bb(end) = Ti;

    TT = AA\bb;
    
    qqp(jj) = hh*(TT(1)-Te);
    
    if jj == 1
        Tis = TT;
    end
end
%% plotto temperatura
plot([.1 .1], [-5 20], 'k-.', 'linewidth', 2)
hold on
plot(xx,Tis,'b-', 'linewidth', 2)
grid on
legend('interfaccia', 'soluzione numerica', 'location', 'southeast')
xlabel('x [m]')
ylabel('T [°C]')
set(gca, 'fontsize', 18)

%% plotto qq
figure()
plot(dkv,qqp,'b-s', 'linewidth', 2)
grid on
legend('perdite', 'location', 'southeast')
ylabel('Q [W/m^2]')
xlabel('k [W/mK]')
set(gca, 'fontsize', 18)
