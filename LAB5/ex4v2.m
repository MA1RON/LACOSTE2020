clc; clear; close all;
%% dati
% geometrici
LL = .5;
bb = .1;
aa = 2e-2;

L0 = .56;

surf = aa*bb;
Awet = 2*(aa+bb);
vol = surf*LL;

% fisici
kk = 350;
Tinf = 4.5;
hh = 500;
II = 7e3;
roel0 = 1.75e-8;

roel = @(z) roel0 * cos(pi*z/L0).^2;
%% ciclo per errori
dzv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2];
err = zeros(size(dzv));
errcond = zeros(size(dzv));
for jjerr = 1:length(dzv)
    %% griglia spaziale
    dz = dzv(jjerr);
    zz = (0:dz:LL/2)'; % dominio simmetrico
    nn = length(zz);

    % calori
    qqel = roel(zz)*LL/surf*II^2/vol;
    qqconv = hh*Tinf*Awet/vol;
    qqq = qqel + qqconv;

    %% costruisco il sistema
    diag_si = ones(nn,1);
    diag_p = -2*ones(nn,1)-hh*Awet/vol*dz^2/kk;
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
    bb = -qqq*dz^2/kk;

    % --- condizioni al contorno ---
    % neumann omogeneo (simmetria)
    AA(1,1) = 1;
    AA(1,2) = -1;
    bb(1) = 0;

    % dirichlet esterno
    AA(end,end-1) = 0;
    AA(end) = 1;
    bb(end) = Tinf;

    %% risolvo il sistema
    TT = AA\bb;
    
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    if jjerr == 1 % soluzione più precisa
        TTref = TT;
        zzref = zz;
        dzref = dz;
    else
        err(jjerr) = norm(TTref(1:round(dz/dzref):end)-TT)/norm(TTref(1:round(dz/dzref):end)-Tinf);
    end
end

%% post production
% andamento temperatura
subplot(1,2,1)
plot(zzref*1e3,TTref,'b-','linewidth',2)
grid on
box on
xlabel('z direction [mm]')
ylabel('temperature [K]')
set(gca,'fontsize',16)

% andamento errore
subplot(1,2,2)
loglog(dzv,err,'b-.*','linewidth',2)
hold on
loglog(dzv,errcond,'r-.*','linewidth',2)
grid on
box on
xlabel('\deltaz [m]')
ylabel('Errore [-]')
set(gca,'fontsize',16)