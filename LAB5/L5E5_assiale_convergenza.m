clc; 
clear; 
close all;

%% dati
% fisici
Tco2 = 300;         % °C
hh = 50;            % W/m^2K
qvol = 8e4;         % W/m^3
kbarra = 350;       % W/mK
kss = 50;           % W/mK

% topologia
d1 = .2;            % m
d2 = .23;           % m
LL = .5;            % m

surfi1 = pi*d1^2/4;         % m^2
surfi2 = pi*(d2^2-d1^2)/4;  % m^2
surfe = pi*(d1 + d2)*LL/2;  % m^2 COSTANTE

vol1 = surfi1*LL;           % m^3
vol2 = surfi2*LL/2;         % m^3
voltot = vol1+vol2;         % m^3 COSTANTE

%% ciclo errori
dxv=[1e-5 5e-5 1e-4 5e-4 1e-3 5e-3];

err = ones(size(dxv));
RR = ones(size(dxv));

for jj = 1:length(dxv)
    %% griglie
    dx = dxv(jj);          
    x1 = (0:dx:LL/2)';    
    x2 = (LL/2+dx:dx:LL)';
    xx = [x1;x2];
    n1 = length(x1);
    n2 = length(x2);
    nn = length(xx);

    keq = (kbarra*(surfi1)+kss*(surfi2))/(surfi1+surfi2);
    kk=kbarra*(xx<=LL/2)+keq*(xx>LL/2);

    %% AA & bb
    diag_i = ones(nn,1);
    diag_p = -2+surfe/voltot*dx^2./-kk*hh;
    diag_s = ones(nn,1);

    AA = spdiags([diag_s diag_p diag_i], -1:1, nn, nn);

    bb = qvol*dx^2./-kk+surfe/voltot*hh*Tco2*dx^2./-kk;

    %% coco & TT
    % robin
    AA(1,1) = -1+hh*dx/-kbarra;
    AA(1,2) = 1;
    bb(1) = hh/-kbarra*Tco2*dx;

    % continuità
    AA(n1,n1-1) = -kbarra;
    AA(n1,n1) = kbarra+keq;
    AA(n1,n1+1) = -keq;
    bb(n1) = 0;

    % robin
    AA(end,end-1) = -1;
    AA(end,end) = 1-hh/-keq*dx;
    bb(end) = hh/-keq*-Tco2*dx;

    TT = AA\bb;
    
    %% errori
    RR(jj) = condest(AA)/(1-condest(AA)*eps)*2*eps;
    
    if jj == 1
        TTref = TT;
        dxref = dx;
        xxref = xx;
    else
        err(jj) = norm(TTref(1:round(dx/dxref):end)-TT)/norm(TTref(1:round(dx/dxref):end)-Tco2);
    end
end

%% plotto errori
figure(1)
loglog(dxv(2:end),err(2:end),'b-o','linewidth',2)
hold on
loglog(dxv,RR,'r-.o','linewidth',2)
hold on
loglog(dxv(3), err(3), 'ko', 'markerfacecolor', 'y','markersize',8)
legend('Troncamento', 'Condizionamento', 'location', 'south')
grid on
box on
title ('Andamento errori')
xlabel('\deltax (m)')
ylabel('Errori (1)')
set(gca,'fontsize',18)