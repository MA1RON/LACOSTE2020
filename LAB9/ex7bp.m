%% LaCoSTe 2020/2021 
%% Lab 9 Ex 7 - Soluzioni, Gianmarco Lorenti

clc;clear all;close all;

set(0,'defaultlinelinewidth',3);
set(0,'defaulttextfontsize',17);
set(0,'defaultaxesfontsize',17);

%% Dati in input

% Barra 
DD = 10e-3; %Diametro della barra (m)
RR = DD/2; %Raggio della barra (m) - asse r

LL = 100e-3; %Lunghezza della barra (m) - asse z

kk = 14; %Conducibilità termica (W/m K)

Tb = 100 + 273.15; %Temperatura della "base" (K)

% Aria (scambio termico convettivo)
Tinf = 25 + 273.15; %Temperatura indisturbata dell'aria (K)
hh_T = @(T) (2.89*(0.6+0.624*(abs(T-Tinf))^(1/6))^2); %Coefficiente di scambio termico (W/m2 K)

% Pareti (scambio termico radiativo)
eps = 0.2; %emissività della superficie della barra (-)
Tsup = 25 + 273.15; %Temperatura della parete (K)
kb=5.670373e-8; %Costante di Boltzmann (W/m2 K4)

qrad_T = @(T) (eps*kb*(Tsup^4 - T.^4)); %Flusso termico radiativo (W/m2)


%% Geometria e discretizzazione del dominio
% Per la simmetria del problema decido di dividere il dominio a metà lungo
% la coordinata radiale, considerando solo un raggio. Dovrò tenerne conto
% andando ad impostare una condizione al contorno di Neumann omogenea
% lungo il bordo sud.
% Creo una griglia rettangolare dove, in prima approssimazione scelgo i
% nodi lungo la coordinata radiale e quella assiale in modo che dr sia più
% grande di dz in modo da cogliere variazioni lungo il raggio, che è molto
% più piccolo della lunghezza.
% Faccio variare la coordinata assiale lungo zvec e quella radiale lungo
% rvec.

nr =51; dr = RR/(nr-1); fprintf('dr: %f [m]\n',dr);
nz = 101; dz = LL/(nz-1); fprintf('dz: %f [m]\n',dz);

ntot = nr*nz; %Numero totale di nodi

rvec = linspace(0,RR,nr); %Considero solo metà dominio, per simmetria
zvec = linspace(0,LL,nz); 

% Creo le matrici che mi serviranno per rappresentare graficamente i
% risultati.
[rmat,zmat] = meshgrid(rvec,zvec);

%% Soluzione del problema: impostazione
% Il problema è non lineare, per cui dovrò risolvere iterativamente,
% impostando un ciclo while. Definisco una tolleranza desiderata
% sull'errore relativo, che definisco come la norma della differenza tra
% TTold e TT [new], che sono due vettori colonna di lunghezza ntot, dove è
% salvata la temperatura in ogni nodo ottenuto dalla "linearizzazione" del
% dominio (dove con linearizzazione intendo il passaggio da una matrice in
% cui ogni nodo è definito da una coppia di indici, (ized,irad), ad un 
% vettore dove ogni nodo è definito da un unico indice, k).

% Per semplificare la scrittura del codice, scelgo di risolvere il problema
% non-lineare sfruttando il metodo dei frozen-coefficient, ovvero vado a
% fissare, ad ogni iterazione, il coefficiente di scambio termico (hh) e il
% calore scambiato per radiazione (qrad), basandomi sulla temperatura
% dell'iterazione precedente. Sempre per semplificare la scrittura, anche
% se ciò richiede di richiamare più volte le funzioni hh_T(T) e qrad_T(T)ad 
% ogni ciclo, decido di andare a calcolare hh e qrad in ogni nodo di 
% interesse.

% Decido infine di tenere tutto ciò che è costante fuori dal ciclo, ovvero
% la matrice dei coefficienti AA e il vettore dei termini noti bb (in
% questo caso fatto di zeri, mancando la generazione interna di calore),
% nonché le modifiche agli array precedenti introdotte dalla condizione al
% contorno di Dirichlet lungo il bordo ovest e quella di Neumann omogenea
% lungo il bordo sud.

% Infine, relativamente alle condizioni al contorno, decido di far
% rientrare tutto il bordo ovest (compresi i vertici nord-ovest e sud-ovest)
% nella condizione di Dirichlet, mentre per i vertici nord-est e sud-est vado 
% a definire condizioni ad-hoc. Riguardo al vertice sud-ovest, imposto la
% condizione al contorno in modo che tenga conto del flusso convettivo e
% radiativo proveniente da est, e del flusso nullo a sud dovuto alla
% simmetria.

% Come ultimo richiamo, nelle equazioni scriverò i vari termini in modo da
% tenere ben distinte le diverse fonti (ad esempio scambio convettivo e
% radiativo), anche se questo dovesse significare costringere il programma
% a fare più operazioni non necessarie.


%% Creo la matrice dei coefficienti
% Attenzione al fatto che ci troviamo in coordinate cilindriche, per cui i
% termini dei coefficienti cambieranno!

AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

for irad = 2:nr - 1
    for ized = 2:nz - 1
        
        r = (irad-1)*dr; %Raggio al nodo (ized,irad)
        
        kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
        kn = kc + nz; %Nodo a nord (ized,irad+1)
        ks = kc - nz; %Nodo a sud (ized,irad-1)
        ke = kc + 1; %Nodo ad est(ized+1,irad)
        kw = kc - 1; %Nodo ad ovest(ized-1,irad)
        
        %fprintf('   %d\n%d %d %d\n   %d\n\n',kn,kw,kc,ke,ks);
        
        AA(kc,kc) = -2*kk*(1/dr^2 + 1/dz^2); %Nodo centrale (ized,irad)
        AA(kc,kw) = kk/dz^2; %Nodo ad ovest(ized-1,irad)
        AA(kc,ke) = kk/dz^2; %Nodo ad est(ized+1,irad)
        AA(kc,kn) = kk*(1 + dr/r)/dr^2; %Nodo a nord (ized,irad+1)
        AA(kc,ks) = kk*(1 - dr/r)/dr^2; %Nodo a sud (ized,irad-1)

    end
end

%% Condizioni al contorno e ciclo while
% Mi muovo in verso antioratio, partendo dal bordo ovest. In ogni caso,
% valgono le osservazioni fatte precedentemente sulle condizioni al
% contorno.

% Bordo ovest (compresi i vertici sud-ovest e nord-ovest): Dirichlet
ized = 1;
for irad = 1 : nr
    
    kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
    AA(kc,kc) = 1;
    bb(kc) = Tb;
    
end

% Bordo sud: Neumann omogenea
irad = 1;
for ized = 2 : nz - 1
    
    kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
    kn = kc + nz; %Nodo a nord (ized,irad+1)
    ke = kc + 1; %Nodo ad est(ized+1,irad)
    kw = kc - 1; %Nodo ad ovest(ized-1,irad)
    
    AA(kc,kc) = -2*(1/dr^2 + 1/dz^2);
    AA(kc,kw) = 1/dz^2; %Nodo ad ovest(ized-1,irad)
    AA(kc,ke) = 1/dz^2; %Nodo ad est(ized+1,irad)
    AA(kc,kn) = 2/dr^2; %Nodo a nord (ized,irad+1)
    
    bb(kc) = 0;
    
end


% Inizializzo le variabili necessarie per il ciclo while
% Per il guess iniziale della temperatura, opto per una temperatura
% uniforme, pari a Tb per tutta la barra. Faccio questa ipotesi in virtù
% della conducibilità relativamente alta della barra e dela piccola
% influenza che mi aspetto dallo scambio termico radiativo e convettivo.
T0 = Tb*ones(size(bb)); %Temperatura di primo tentativo

toll = 1e-4; %Tolleranza sull'errore relativo (-)
err = toll + 1; %Errore relativo (-) 
cont = 0; %Contatore delle iterazioni

TTold = Tb*ones(ntot,1); %Temperatura dell'iterazione precedente

while err > toll
    
    cont = cont +1;
    
    % Continuo con l'assegnazione delle condizioni al contorno, per quelle
    % condizioni che rendono il problema non-lineare e pertanto hanno
    % termini che variano con la temperatura.
    % Attenzione, decido di assegnare tutta la condizione al contorno
    % all'interno del ciclo ai fini di rendere più chiaro quanto fatto, ma
    % potrei anche assegnare solo le parti che effettivamente variano
    % (AA(kc,kc) e bb(kc)).
    
    
    % Vertice sud-est: Robin + Neumann ad ovest / Neumann omogenea a sud
    irad = 1;
    ized = nz;
    
    kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
    kn = kc + nz; %Nodo a nord (ized,irad+1)
    kw = kc - 1; %Nodo ad ovest(ized-1,irad)
    
    % Attenzione, decido di calcolare h e qrad per ogni nodo interessato,
    % usando le funzioni definite in precedenza. Questi valori sono
    % calcolati usando la temperatura dell'iterazione precedente
    
    Tkc = TTold(kc);
    hh = hh_T(Tkc);
    qrad = qrad_T(Tkc);
    
    AA(kc,kc) = -(1/dr^2 + 1/dz^2 + hh/(kk*dz) + 0);
    AA(kc,kw) = 1/dz^2; %Nodo ad ovest(ized-1,irad)
    AA(kc,kn) = 1/dr^2; %Nodo a sud (ized,irad-1)
    
    bb(kc) = -(hh/(kk*dz) + 0)*Tinf -(1/(kk*dz) + 0)*qrad;
    
    
    % Bordo est: Robin + Neumann
    ized = nz;
    for irad = 2 : nr - 1
        
        kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
        kn = kc + nz; %Nodo a nord (ized,irad+1)
        ks = kc - nz; %Nodo a sud (ized,irad-1)
        kw = kc - 1; %Nodo ad ovest(ized-1,irad)
        
        Tkc = TTold(kc);
        hh = hh_T(Tkc);
        qrad = qrad_T(Tkc);
        
        AA(kc,kc) = -2*(1/dr^2 + 1/dz^2 + hh/(kk*dz));
        AA(kc,kw) = 2/dz^2; %Nodo ad ovest(ized-1,irad)
        AA(kc,kn) = 1/dr^2; %Nodo a nord (ized,irad+1)
        AA(kc,ks) = 1/dr^2; %Nodo a sud (ized,irad-1)
        
        bb(kc) = -2*hh/(kk*dz)*Tinf -2*qrad/(kk*dz);
        
    end
    
    
    % Vertice nord-est: Robin + Neumann a ovest / Robin + Neumann a nord
    irad = nr;
    ized = nz;
    
    kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
    ks = kc - nz; %Nodo a sud (ized,irad-1)
    kw = kc -1; %Nodo ad ovest(ized-1,irad)
    
    Tkc = TTold(kc);
    hh = hh_T(Tkc);
    qrad = qrad_T(Tkc);
    
    AA(kc,kc) = -(1/dr^2 + 1/dz^2 + hh/(kk*dr) + hh/(kk*dz));
    AA(kc,kw) = 1/dz^2; %Nodo ad ovest(ized-1,irad)
    AA(kc,ks) = 1/dr^2; %Nodo a sud (ized,irad-1)
    
    bb(kc) = -(hh/(kk*dr) + hh/(kk*dz))*Tinf -(1/(kk*dr) + 1/(kk*dz))*qrad;
   
        
    % Bordo nord: Robin + Neumann
    irad = nr;
    for ized = 2 : nz - 1
        
        kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
        ks = kc - nz; %Nodo a sud (ized,irad-1)
        ke = kc + 1; %Nodo ad est(ized+1,irad)
        kw = kc - 1; %Nodo ad ovest(ized-1,irad)
        
        Tkc = TTold(kc);
        hh = hh_T(Tkc);
        qrad = qrad_T(Tkc);
        
        AA(kc,kc) = -2*(1/dr^2 + 1/dz^2 + hh/(kk*dr));
        AA(kc,kw) = 1/dz^2; %Nodo ad ovest(ized-1,irad)
        AA(kc,ke) = 1/dz^2; %Nodo ad est(ized+1,irad)
        AA(kc,ks) = 2/dr^2; %Nodo a sud (ized,irad-1)
        
        bb(kc) = -2*hh/(kk*dr)*Tinf -2*qrad/(kk*dr);
        
    end


    TT = AA\bb; %Temperatura dell'iterazione corrente (K), soluzione numerica
    
    err = norm(TT-TTold)/norm(TT-Tinf); %Calcolo l'errore relativo
    
    TTold = TT; %Aggiorno la temperatura dell'iterazione precedente


    % Salvo in un array la temperatura media della punta (bordo est).
    % Potrei farlo anche solo alla fine del while, quando si raggiunge la
    % convergenza, ma decido piuttosto di salvarne il valore ad ogni
    % iterazione in modo da vederne l'evoluzione
    Tm = 0;
    contTm = 0;
    
    ized = nz;
    for irad = 1 : nr 
        
        kc = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)
        Tm = Tm + TT(kc);
        contTm = contTm + 1;
        
    end
    
    TTm(cont) = Tm/contTm;

end


%% Post-process

% Distribuzione di temperatura 2D
figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);
f1 = figure(1);

TTplot = reshape(TT,nz,nr);
surf(zmat*1e3,rmat*1e3,TTplot-273.15);

xlabel('z (mm)');
ylabel('r (mm)');
zlabel('T (°C)');
colorbar;

print -f1 -depsc Lab9Ex7_figure1


% Distribuzione di temperatura 1D
figure('units','normalized','outerposition',[0 0 0.5 0.5]);
f2=figure(2);

hold on;
grid on;
box on;

izplot = round(linspace(1,nz,5));

for ii = 1:length(izplot)
    
    ized = izplot(ii);
    
    for irad = 1 : nr      
    kkc(irad) = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)      
    end
    
    Trplot = TT(kkc);
    
    plot(rvec*1e3,Trplot - 273.15,'displayname',sprintf('z: %.1f mm',zvec(ized)*1e3));
    
end

xlabel('r (mm)');
ylabel('T (°C)');

ylim([(Tsup - 273.15), 1.2*(Tb - 273.15)]);

legend('-dynamiclegend','location','eastoutside','orientation','vertical');
title('Temperatura lungo il raggio, a diversi z');


% Distribuzione di temperatura 1D
zplot = input('\nInserisci un valore della coordinata assiale: ');
while zplot < 0 || zplot > LL*1e3 
    zplot =  input(sprintf('Il valore deve essere compreso tra %f e %f: ',0,LL*1e3)); 
end
[~ , izplot] = min( abs( zvec*1e3 - zplot ) );
zplot = zvec(izplot)*1e3;

ized = izplot;
for irad = 1 : nr      
    kkc(irad) = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)      
end

Trplot = TT(kkc);


rplot = input('\nInserisci un valore della coordinata radiale: ');
while rplot < 0 || rplot > RR*1e3
    rplot =  input(sprintf('Il valore deve essere compreso tra %f e %f: ',0,RR*1e3)); 
end
[~ , irplot] = min( abs( rvec*1e3 - rplot ) );
rplot = rvec(irplot)*1e3;

irad = irplot;
for ized = 1 : nz      
    kkc(ized) = (irad - 1)*nz + ized; %Nodo centrale (ized,irad)      
end

Tzplot = TT(kkc);

figure('units','normalized','outerposition',[0.5 0 0.5 1]);
f3=figure(3);

subplot(2,1,1)

hold on;
grid on;
box on;

plot(rvec*1e3,Trplot - 273.15,'b');
xlabel('r (mm)');
ylabel('T (°C)');

title(sprintf('(a): Temperatura lungo il raggio, a z pari a %.1f mm',zplot));

subplot(2,1,2)

hold on;
grid on;
box on;

plot(zvec*1e3,Tzplot - 273.15,'r');
xlabel('z (mm)');
ylabel('T (°C)');

title(sprintf('(a): Temperatura parallelamente all''asse, a r pari a %.1f mm',rplot));

print -f3 -depsc Lab9Ex7_figure3

fprintf('\nLa temperatura media della punta è pari a: %.1f °C\n',TTm(end) - 273.15);

