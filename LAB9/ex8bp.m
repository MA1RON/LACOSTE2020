%% LaCoSTe 2020/2021 
%% Lab 9 Ex 7 - Soluzioni, Gianmarco Lorenti

clc;clear all;close all;

set(0,'defaultlinelinewidth',3);
set(0,'defaulttextfontsize',17);
set(0,'defaultaxesfontsize',17);

% Una piastra (k = 10 [W/m/K]) è rinforzata da una serie di nervature 
% longitudinali a sezione rettangolare di lunghezza L = 8 [mm] e spessore
% w = 4 [mm]. La base del piatto è mantenuta alla temperatura costante 
% Tb = 45°C, mentre le superfici delle nervature sono esposte all’aria a 
% temperatura Tinf = 25°C e coefficiente di scambio termico h = 600 [W/m2/K].
 
% a. Usando il metodo delle differenze finite calcola la distribuzione di 
%    temperatura e il calore scambiato alla base della nervatura.
% b. Valuta l’accuratezza della soluzione con uno studio di grid 
%    independence

%% Dati in input

% Nervature
LL = 8e-3; %Lunghezza (m)
ww = 4e-3; %Spessore (m)

kk = 10; %Conducibilità termica (W/m K)

Tb = 45 + 273.15; %Temperatura alla base (K)

% Scambio termico con l'aria
hh = 600; %Coefficiente di scambio termico co l'aria (W/m2 K)
Tinf = 25 + 273.15; %Temperatura indisturbata dell'aria (K)

%% Geometria e discretizzazione del dominio
% Per la simmetria del problema decido di dividere il dominio a metà lungo
% la coordinata dello spessore (y), considerandone solo metà. Dovrò tenerne 
% conto andando ad impostare una condizione al contorno di Neumann omogenea
% lungo il bordo sud.
% Creo una griglia rettangolare dove, in prima approssimazione scelgo i
% nodi lungo la coordinata y e quella x in modo che dx e dy abbiano lo
% stesso valore, visto che lo spessore e la lunghezza sono paragonabili.
% Faccio variare la coordinata x lungo xvec e quella y lungo yvec.

nx = 21; dx = LL/(nx-1); fprintf('dx: %f [m]\n',dx);
ny = 6; dy = ww/2/(ny-1); fprintf('dy: %f [m]\n',dy);

ntot = ny*nx; %Numero totale di nodi

xvec = linspace(0,LL,nx); 
yvec = linspace(0,ww/2,ny); %Considero solo metà dominio, per simmetria

% Creo le matrici che mi serviranno per rappresentare graficamente i
% risultati.
[xmat,ymat] = meshgrid(xvec,yvec);


%% Soluzione del problema: impostazione
% Il problema è lineare, quindi può essere risolto in modo relativamente
% facile. Devo scrivere una matrice AA dei coefficienti e un vettore bb dei
% termini noti che tengano conto dell'equazione di conduzione del calore
% (in due dimensioni) in ogni nodo e delle condizioni al contorno lungo i
% bordi.

% Relativamente a queste ultime, decido di far rientrare tutto il bordo
% ovest (compresi i vertici nord-ovest e sud-ovest)nella condizione di 
% Dirichlet, mentre per i vertici nord-est e sud-est vado a definire 
% condizioni ad-hoc. Riguardo al vertice sud-ovest, imposto la condizione 
% al contorno in modo che tenga conto del flusso convettivo e radiativo
% proveniente da est, e del flusso nullo a sud dovuto alla simmetria.

% Come ultimo richiamo, nelle equazioni scriverò i vari termini in modo da
% tenere ben distinte le diverse fonti (ad esempio scambio convettivo e
% radiativo), anche se questo dovesse significare costringere il programma
% a fare più operazioni non necessarie.


%% Creo la matrice dei coefficienti
% Attenzione al fatto che ci troviamo in coordinate cartesiane, per cui i
% termini dei coefficienti cambieranno!

AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

for iy = 2:ny- 1
    for ix = 2:nx - 1
        
        kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
        kw = kc - 1; %Nodo ad ovest(ix-1,iy)
        ks = kc - nx; %Nodo a sud (ix,iy-1)
        ke = kc + 1; %Nodo ad est(ix+1,iy)
        kn = kc + nx; %Nodo a nord (ix,iy+1)
        
        AA(kc,kc) = -2*kk*(1/dx^2 + 1/dy^2); %Nodo centrale (ix,iy)
        AA(kc,kw) = kk/dx^2; %Nodo ad ovest(ix-1,iy)
        AA(kc,ks) = kk/dy^2; %%Nodo a sud (ix,iy-1)
        AA(kc,ke) = kk/dx^2; %Nodo ad est(ix+1,iy)
        AA(kc,kn) = kk/dy^2; %Nodo a nord (ix,iy+1)
        
    end
end

%% Condizioni al contorno 
% Mi muovo in verso antioratio, partendo dal bordo ovest. In ogni caso,
% valgono le osservazioni fatte precedentemente sulle condizioni al
% contorno.


% Bordo ovest (compresi i vertici sud-ovest e nord-ovest): Dirichlet
ix = 1;
for iy = 1 : ny
    
    kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
    AA(kc,kc) = 1;
    bb(kc) = Tb;
    
end


% Bordo sud: Neumann omogenea
iy = 1;
for ix = 2 : nx - 1
    
    kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
    kw = kc - 1; %Nodo ad ovest(ix-1,iy)
    ke = kc + 1; %Nodo ad est(ix+1,iy)
    kn = kc + nx; %Nodo a nord (ix,iy+1)
   
    AA(kc,kc) = -2*(1/dx^2 + 1/dy^2);
    AA(kc,kw) = 1/dx^2; %Nodo ad ovest(ix-1,iy)
    AA(kc,ke) = 1/dx^2; %Nodo ad est(ix+1,iy)
    AA(kc,kn) = 2/dy^2; %Nodo a nord (ix,iy+1)
    
    bb(kc) = 0;
    
end


% Vertice sud-est: Robin ad ovest / Neumann omogenea a sud
iy = 1;
ix = nx;

kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
kw = kc - 1; %Nodo ad ovest(ix-1,iy)
kn = kc + nx; %Nodo a nord (ix,iy+1)

AA(kc,kc) = -(1/dx^2 + 1/dy^2 + hh/(kk*dx) + 0);
AA(kc,kw) = 1/dx^2; %Nodo ad ovest(ix-1,iy)
AA(kc,kn) = 1/dy^2; %Nodo a nord (ix,iy-1)

bb(kc) = -(hh/(kk*dx) + 0)*Tinf;


% Bordo est: Robin 
ix = nx;
for iy = 2 : ny - 1
    
     kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
     kw = kc - 1; %Nodo ad ovest(ix-1,iy)
     ks = kc - nx; %Nodo a sud (ix,iy-1)
     kn = kc + nx; %Nodo a nord (ix,iy+1)
    
     AA(kc,kc) = -2*(1/dx^2 + 1/dy^2 + hh/(kk*dx));
     AA(kc,kw) = 2/dx^2; %Nodo ad ovest(ix-1,y)
     AA(kc,ks) = 1/dy^2; %Nodo a sud (ix,iy-1)
     AA(kc,kn) = 1/dy^2; %Nodo a nord (ix,iy+1)
    
     bb(kc) = -2*hh/(kk*dx)*Tinf;
    
end


% Vertice nord-est: Robin a ovest / Robin a nord
iy = ny;
ix = nx;

kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
kw = kc -1; %Nodo ad ovest(ix-1,iy)
ks = kc - nx; %Nodo a sud (ix,iy-1)

AA(kc,kc) = -(1/dx^2 + 1/dy^2 + hh/(kk*dx) + hh/(kk*dy));
AA(kc,kw) = 1/dx^2; %Nodo ad ovest(ix-1,iy)
AA(kc,ks) = 1/dy^2; %Nodo a sud (ix,iy-1)

bb(kc) = -(hh/(kk*dx) + hh/(kk*dy))*Tinf;


% Bordo nord: Robin
iy = ny;
for ix = 2 : nx - 1
    
    kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
    kw = kc - 1; %Nodo ad ovest(ix-1,iy)
    ks = kc - nx; %Nodo a sud (ix,iy-1)
    ke = kc + 1; %Nodo ad est(ix+1,iy)

    AA(kc,kc) = -2*(1/dx^2 + 1/dy^2 + hh/(kk*dy));
    AA(kc,kw) = 1/dx^2; %Nodo ad ovest(ix-1,iy)
    AA(kc,ks) = 2/dy^2; %Nodo a sud (ix,iy-1)
    AA(kc,ke) = 1/dx^2; %Nodo ad est(ix+1,iy)
    
    bb(kc) = -2*hh/(kk*dy)*Tinf;
    
end

TT = AA\bb; %Temperatura dell'iterazione corrente (K), soluzione numerica


ix = 1;

iy = 1;
kc = (iy - 1)*nx + ix; ke = kc + 1; kn = kc + nx; %Nodi centrale, est, nord
qbase = 0 + kk*(dy*(TT(ke)-TT(kc))/dx + dx/2*(TT(kn)-TT(kc))/dy); %Calore scambiato al nodo sud (W/m)

for iy = 2 : ny - 1
    
    kc = (iy - 1)*nx + ix; %Nodo centrale (ix,iy)
    ks = kc - nx; %Nodo a sud(ix,iy-1)
    ke = kc + 1; %Nodo ad est(ix+1,iy)
    kn = kc + nx; %Nodo a nord (ix,iy+1)
    
    qbase = qbase + kk*(dx/2*(TT(ks)-TT(kc))/dy + dy*(TT(ke)-TT(kc))/dx +...
        dx/2*(TT(kn)-TT(kc))/dy); %Calore scambiato ad ogni nodo (W/m)
    
end

iy = ny;
kc = (iy - 1)*nx + ix; ks = kc - nx; ke = kc + 1; %Nodi centrale, sud, est
qbase = qbase + kk*(dx/2*(TT(ks)-TT(kc))/dy + dy*(TT(ke)-TT(kc))/dx); %Calore scambiato al nodo sud (W/m)


%% Post-process

% Distribuzione di temperatura 2D
figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);
f1 = figure(1);

TTplot = reshape(TT,nx,ny);
surf(xmat*1e3,ymat*1e3,TTplot'-273.15);

xlabel('z (mm)');
ylabel('r (mm)');
zlabel('T (°C)');
colorbar;

print -f1 -depsc Lab9Ex8_figure1
