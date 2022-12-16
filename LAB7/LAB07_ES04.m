%% LaCoSTe 2020/2021 
%% Lab 7 Ex 4 - Soluzioni, Gianmarco Lorenti

clc;clear all;close all;

set(0,'defaultlinelinewidth',3);
set(0,'defaulttextfontsize',17);
set(0,'defaultaxesfontsize',17);

%% Input

% Delle sfere di acciaio di diametro 12 mm (conducibilità = 40 W/m/K, 
% densità = 7800 kg/m3, e calore specifico = 600 J/kg/K) vengono temprate 
% scaldandole molto rapidamente a 1150 K e poi raffreddandole lentamente
% fino alla temperatura finale di 400K in un ambiente con aria, la cui 
% temperatura Taria  aumenta nel tempo:  Taria = 325 K + 0.1875 K/s × t, 
% dove t è il tempo dall’inizio del processo di raffreddamento.
% Calcolare numericamente con lo schema di Eulero implicito l’evoluzione 
% nel tempo della temperatura delle sferette, supponendo che il coefficiente
% di scambio termico sia pari a 20 W/m2K e che tutta la superficie esterna 
% delle sfere sia in contatto con l’aria. Dopo aver fatto il grafico 
% dell’evoluzione nel tempo della temperatura delle sfere, determinare dopo 
% quanto tempo essa diventa uguale a quella dell’aria. Verificare il 
% risultato ottenuto utilizzando lo schema di Eulero in avanti.

% Sfere acciaio

phi = 12e-3; %Diametro sfere (mm)
kk = 40; %Conducibilità termica acciacio (W/m K)
rho = 7800; %Densità acciacio (kg/m3)
cc = 600; %Calore specificio acciaio (J/kg K)

% Processo di tempra

Tin = 1150; %Temperatura iniziale delle sfere (K)
Tf = 400; %Temperatura che dovrebbero raggiungere le sfere al termine del processo (K)

Ta = @(t) 325 + 0.1875*t; %Temperatura dell'aria in funzione del tempo (K), t(s)
hh = 20; %Coefficiente di scambio termico acciaio-aria (W/m2 K)

% Numero di Biot

Bi = hh*phi/kk; %Analisi monodimensionale quando Bi << 0.1

%% Analisi 0D - BE e FE

T0 = Tin; %"Distribuzione" iniziale di temperatura
Vol = 4/3*pi*(phi/2)^3; %Volume delle sfere (m3)
Sup = 4*pi*(phi/2)^2; %Superficie di scambio termico (m2)


% Backward Euler (metodo implicito)

Told_be = T0; %Inizializzazione temperatura sfere
t_be = 0; %Istante di tempo iniziale (s)
dt_be = 10; %Timestep (s)

count_be = 1; %Contatore iterazioni
time_be(count_be) = t_be; %Vettore dei tempi 
Tt_be(count_be) = Told_be; %Vettore della temperatura delle sfere nel tempo
Tat_be(count_be) = Ta(t_be); %Vettore della temperatura dell'aria nel tempo 

% Si avvia un ciclo while per il calcolo della temperatura delle sfere (in
% funzione di quella dell'aria) ad ogni istante di tempo t + dt. Il ciclo si
% arresta quando la temperatura delle sfere raggiunge quella dell'aria. E'
% necessario prestare attenzione a questo passaggio perché la temperatura
% dell'aria è data in funzione del tempo soltanto e quindi non è
% influenzata dalla temperatura delle sfere. In realtà questa è
% un'approssimazione poiché, per il secondo principio, la temperatura
% dell'aria non può superare quella delle sfere.
% Si aggiunge un numero massimo di iterazioni per sicurezza.
% Viene utilizzato un metodo esplicito (Forward Euler).

while (Told_be > Tf && Told_be>Ta(t_be)) && count_be < 1e4
    
    count_be = count_be + 1; %Aggiornamento contatore iterazioni
    t_be = t_be + dt_be; %Aggiornamento tempo
    time_be(count_be) = t_be;
    Tat_be(count_be) = Ta(t_be); %Calcolo temperatura dell'aria
    
    Tnew_be = (Told_be + (hh*Sup*dt_be)/(rho*cc*Vol)*Ta(t_be))/(1+(hh*Sup*dt_be)/(rho*cc*Vol));
    
    Tt_be(count_be) = Tnew_be;
    Told_be = Tnew_be; %Aggiornamento di Told ad ogni iterazione!!
    
end
  

% Forward Euler (metodo esplicito)

Told_fe = T0; %Inizializzazione temperatura sfere
t_fe = 0; %Istante di tempo iniziale (s)
dt_fe = 10; %Timestep (s)

count_fe = 1; %Contatore iterazioni
time_fe(count_fe) = t_fe; %Vettore dei tempi
Tt_fe(count_fe) = Told_fe; %Vettore della temperatura delle sfere nel tempo
Tat_fe(count_fe) = Ta(t_fe); %Vettore della temperatura dell'aria nel tempo

% Si avvia un ciclo while per il calcolo della temperatura delle sfere (in
% funzione di quella del gas) ad ogni istante di tempo t + dt. Il ciclo si
% arresta quando la temperatura delle sfere raggiunge quella dell'aria. E'
% necessario prestare attenzione a questo passaggio perché la temperatura
% dell'aria è data in funzione del tempo soltanto e quindi non è
% influenzata dalla temperatura delle sfere. In realtà questa è
% un'approssimazione poiché, per il secondo principio, la temperatura
% dell'aria non può superare quella delle sfere.
% Si aggiunge un numero massima di iterazioni per sicurezza.
% Viene utilizzato un metodo esplicito (Forward Euler).

while (Told_fe > Tf && Told_fe>Ta(t_fe)) && count_fe < 1e4
    
    count_fe = count_fe + 1; %Aggiornamento contatore iterazioni
    t_fe = t_fe + dt_fe; %Aggiornamento tempo
    time_fe(count_fe) = t_fe;
    Tat_fe(count_fe) = Ta(t_fe); %Calcolo temperatura dell'aria
    
    Tnew_fe = Told_fe + (hh*Sup)/(rho*cc*Vol)*dt_fe*(Ta(t_fe) - Told_fe);
    
    Tt_fe(count_fe) = Tnew_fe;
    Told_fe = Tnew_fe; %Aggiornamento di Told ad ogni iterazione!!
    
end


%% Post-process

% Distribuzione di temperatura
figure('units','normalized','outerposition',[0.25 0.5 0.5 0.5]);
f1=figure(1);

hold on;
grid on;
box on;

plot(time_be,Tt_be,'b');
plot(time_fe,Tt_fe,'r-.');
plot(time_be,Tat_be,'k','linewidth',2);
plot(time_be(end),Tt_be(end),'rs','markerfacecolor','w','markersize',8);

xlabel('t (s)');
ylabel('T (K)');

xlim([time_be(1),time_be(end)]);

legend('BE','FE','Aria','orientation','horizontal','location','northeast');

% title('Evoluzione della temperatura media della sfera nel tempo');

text(time_be(end)-10,Tt_be(end)+30,sprintf('%d s',time_be(end)),...
    'horizontalalignment','right','verticalalignment','bottom','color','r','edgecolor','k','linestyle','-');

print -f1 -depsc Lab7Ex4_figure1

if Tt_be(end) > Tf
    fprintf('\nIl tempo necessario perché la temperatura delle sfere eguagli quella dell''aria è pari a %d s (BE).\n',...
        time_be(end));
    fprintf('\nTuttavia il raffreddamento non è sufficiente poiché si raggiunge una temperatura di %.1f K invece che %.1f K.\n',...
        Tt_be(end),Tf);
else
    fprintf('\nIl tempo necessario perché la temperatura delle sfere raggiunga la temperatura finale è pari a %d s (BE).\n',...
        time_be(end));
     fprintf('\nTuttavia il raffreddamento non è completato poiché si l''aria è ancora ad una tmperatura di %.1f K.\n',...
        Tat_be(end));
end
    
