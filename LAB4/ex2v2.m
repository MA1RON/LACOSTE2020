clc; clear; close all;
%% dati
% geometrici
DD = .01;
LL = 4;

surf = pi*DD^2/4;
Awet = pi*DD*LL;
vol = surf*LL;

% fisici
II = 1e3;
Tb = 77;
hh = 500;
kk = 350;
roel = 1.75e-8;

%% valuto errore
dxv = [1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2];
err = zeros(size(dxv));
errcond = zeros(size(dxv));
for jjerr = 1:length(dxv)
    %% griglia
    dx = dxv(jjerr);
    xx = (0:dx:LL)';
    nn = length(xx);

    %% creo il problema
    diag_si = ones(nn,1);
    diag_p = -2*ones(nn,1) - hh*Awet/vol/kk*dx^2;
    AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
    bb = -(roel*LL/surf*II^2/vol + hh*Tb*Awet/vol)/kk*ones(nn,1)*dx^2;

    % condizioni al contorno
    AA(1,1) = 1; AA(1,2) = 0; bb(1) = Tb;
    AA(end,end-1) = 0; AA(end,end) = 1; bb(end) = Tb;

    %% risolvo il problema
    TT = AA\bb;
    
    %% valuto errore
    if jjerr == 1 % solluzione più precisa
        dxref = dx;
        xxref = xx;
        TTref = TT;
    else
        err(jjerr) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref(1:round(dx/dxref):end)-Tb);
    end
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
end

%% post production
% temperatura
figure
plot(xx,TT,'b-','linewidth',2)
hold on
plot([xx(1) xx(end)],[Tb Tb],'k-.','linewidth',2)
legend('Cilindro','Fluido & bordi','location','northeastoutside')
title('Stazionario 1D')
grid on
box on
xlabel('Lunghezza [m]')
ylabel('Temperatura [K]')
set(gca,'fontsize',18)

% errore
figure
loglog(dxv,err,'b-.*','linewidth',2)
hold on
loglog(dxv,errcond,'r-.*','linewidth',2)
legend('Errore di discretizzazione','Upper bound di condizionamento','location','southwest')
title('Errore stazionario 1D')
grid on
box on
xlabel('\deltax [m]')
ylabel('Errore [-]')
set(gca,'fontsize',18)

%% conservazione dell'energia
qvol = roel*II^2*LL/surf/vol; % W/m^3
Eint = qvol*vol; % W
Eout = abs(hh*(Tb-max(TT)))*Awet; % W