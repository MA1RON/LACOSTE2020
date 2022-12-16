clc; clear; close all;

%% dati
kcu = 150; %W/mK
cpcu = 350; %J/kgK
rocu = 8900; %kg/m^3
roel = 2e-8; %hom*m
II = 2e4; %A
TN = 77; %K
hi = 1.2e3; %W/m^2K
kis = .5; %W/mK
cpis = 1.8e3; %J/kgK
rois = 2.5e3; %kg/m^3
Ta = 25+273.15; %K
he = 15; %W/m^2K
ri = 15e-3; %m
rcu = 25e-3; %m
re = 26e-3; %m

qqq=II^2*roel/(pi*(rcu^2-ri^2))^2; %W/m^3

%% griglie
drv = sort([logspace(-7,-3,5) 2*logspace(-7,-3,5) 5*logspace(-7,-3,5)]);
nnv = ones(size(drv));
timev = ones(size(drv));
for jj = 1:length(drv)
    dr = drv(jj);
    rrcu = (ri:dr:rcu)';
    nncu = length(rrcu);
    rris = (rcu+dr:dr:re)';
    nnis = length(rris);

    rr = [rrcu; rris];
    nn = nncu+nnis;
    nnv(jj) = nn;

    %% AA & bb
    diag_s = [dr+2*rrcu; dr+2*rris];
    diag_p = [-4*rrcu; -4*rris];
    diag_i = [-dr+2*rrcu; -dr+2*rris];

    AA = spdiags([diag_s diag_p diag_i],-1:1, nn,nn);

    bb = [qqq/-kcu*2*dr^2*rrcu;zeros(nnis,1)];

    %% coco & TT
    % robin interna
    AA(1,1) = 1+hi*dr/kcu;
    AA(1,2) = -1;
    bb(1) = hi*TN/kcu*dr;

    % continuità di flusso
    AA(nncu,nncu-1) = -kcu*dr;
    AA(nncu,nncu) = kcu*dr+kis*dr;
    AA(nncu,nncu+1) = -kis*dr;
    bb(nncu) = 0;

    % robin esterna
    AA(end,end-1) = -1;
    AA(end,end) = 1+he*dr/kis;
    bb(end) = he*Ta*dr/kis;
    
    tic
    
    TT = AA\bb;
    
    timev(jj) = toc;
    
    if jj == 1
        rref = rr;
        Tref = TT;
    end
end
%% plotto T
plot(rref, Tref,'b-','linewidth', 2);
hold on
plot([25e-3 25e-3], [135 144], 'k-.', 'linewidth', 2);
legend('soluzione numerica', 'interfaccia', 'location', 'south')
xlabel('x [m]')
ylabel('T [K]')
grid on

close
%% plotto time
figure()
semilogy(nnv, timev, 'b-s', 'linewidth', 2);
legend('tempo', 'location', 'south')
xlabel('numero nodi [1]')
ylabel('tempo [s]')
grid on