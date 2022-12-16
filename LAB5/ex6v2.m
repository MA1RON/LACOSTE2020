clc; clear; close all;
%% dati
% dominio
% r -> hhint | 0 - ss - n1 - cu - n2 - is - end | hhext
% geometrici
r1 = 10e-3;
r2 = 14e-3;
r3 = 17e-3;
r4 = 18.5e-3;

LL = 1; % solo per i calcoli, si semplifica

surf1 = pi*r1^2;
surf2 = pi*r2^2;
surf3 = pi*r3^2;
surf4 = pi*r4^2;
surftot = surf4-surf1;

% fisici
TTint = 25;
hhint = 2e3;
TText = 15;
hhext = 5;

IInom = 16e3;

kkss = 1.5;
kkcu = 30;
kkis = .21;

roelss = 5e-6;
roelcu = 1.8e-8;
roelis = 1e-2;

RRss = roelss*LL/(surf2-surf1);
RRcu = roelcu*LL/(surf3-surf2);
RRis = roelis*LL/(surf4-surf3);

IIss = IInom*(RRcu*RRis)/(RRss*RRcu+RRcu*RRis+RRis*RRss);
IIcu = IInom*(RRss*RRis)/(RRss*RRcu+RRcu*RRis+RRis*RRss);
IIis = IInom*(RRss*RRcu)/(RRss*RRcu+RRcu*RRis+RRis*RRss);

drv = [1e-5 2e-5 5e-5 1e-4 2e-4 5e-5 1e-3];
err = zeros(size(drv));
errcond = zeros(size(drv));
for jjerr = 1:length(drv)
    %% griglia radiale
    dr = drv(jjerr);
    rr = (r1:dr:r4)';
    nn = length(rr);
    n1 = length(r1:dr:r2);
    n2 = length(r1:dr:r3);

    % variabili in funzione del raggio
    kk = kkss*(rr<=r2) + kkcu*(rr>r2 & rr<=r3) + kkis*(rr>r3 & rr<=r4);

    surfi = (surf2-surf1)*(rr<=r2) + (surf3-surf2)*(rr>r2 & rr<=r3) + (surf4-surf3)*(rr>r3 & rr<=r4);
    roel = roelss*(rr<=r2) + roelcu*(rr>r2 & rr<=r3) + roelis*(rr>r3 & rr<=r4);
    II = IIss*(rr<=r2) + IIcu*(rr>r2 & rr<=r3) + IIis*(rr>r3 & rr<=r4);

    qqq = II.^2.*roel./surfi.^2;

    %% sistema
    diag_s = .5+rr/dr;
    diag_p = -2*rr/dr;
    diag_i = -.5+rr/dr;
    AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1,nn,nn);
    bb = -qqq./kk*dr.*rr;

    % --- condizioni al contorno ---
    % robin interno
    AA(1,1) = -1-hhint/kkss*dr;
    AA(1,2) = 1;
    bb(1) = -hhint/kkss*TTint*dr;

    % contatto ss - cu
    AA(n1,n1-1) = -1;
    AA(n1,n1) = 1+kkcu/kkss;
    AA(n1,n1+1) = -kkcu/kkss;
    bb(n1) = 0;

    % contatto cu - is
    AA(n2,n2-1) = -1;
    AA(n2,n2) = 1+kkis/kkcu;
    AA(n2,n2+1) = -kkis/kkcu;
    bb(n2) = 0;

    % robin esterno
    AA(end,end-1) = 1;
    AA(end,end) = -1-hhext/kkis*dr;
    bb(end) = -hhext*TText/kkis*dr;

    % --- risolvo il sistema ---
    TT = AA\bb;
    
    %% valuto errori
    errcond(jjerr) = condest(AA)/(1-eps*condest(AA))*2*eps;
    if jjerr == 1
        TTref = TT;
        drref = dr;
        rrref = rr;
    else
        err(jjerr) = norm(TTref(1:round(dr/drref):end)-TT)/norm(TTref-TText);
    end
end

%% post production
figure('position', [100 100 1200 700])

subplot(1,2,1)
plot(rrref*1e3,TTref,'linewidth',2)
grid on
box on
title('Distribuzione di temperatura')
xlabel('Raggio (mm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)

subplot(1,2,2)
loglog(drv,errcond,'r-.*','linewidth',2)
hold on
loglog(drv,err,'b-.*','linewidth',2)
legend('Upper bound troncamento','Errore di discretizzazione','location','southwest')
grid on
box on
title('Andamento errore')
xlabel('\Deltar (m)')
ylabel('Errore (-)')
set(gca,'fontsize',22)