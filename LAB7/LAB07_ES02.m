%% LaCoSTe 2020/2021 
%% LAB07_ES02 - Soluzioni, Matteo Nicoli
clear all
close all
clc

%% Dati di input
% Proprietà titanio
d=0.1;
L=0.2;
V_Ti=pi/4*d^2*L;
k_Ti=21.9;
rho_Ti=4500;
cp_Ti=520;
T(1,1)=600+273.15;

% Proprietà olio
V_oil=5e-3;
k_oil=0.1;
rho_oil=900;
cp_oil=7200;
T(2,1)=15+273.15;

% Scambio termico
h=100;
S=pi*d*L+2*pi/4*d^2;

%% Soluzione numerica (proprietà costanti)
tol=1e-4;
err=1;
i=1;

% Numero di Biot e dt_lim
L=V_Ti/S;
Bi=h*L/k_Ti;
tau=rho_Ti*V_Ti*cp_Ti/h/S;

% Imposto il dt e l'istante iniziale
dt=10;
t=0;

% Per l'uscita dal ciclo while, ho scelto di verificare quando entrambi gli
% errori relativi (quello associato alla temperatura dell'olio e quello
% associato alla temperatura del titanio) sono minori della tollerenza.
coeff_oil=V_oil*rho_oil*cp_oil/h/S;
coeff_Ti=V_Ti*rho_Ti*cp_Ti/h/S;
A=[1-dt/coeff_Ti,dt/coeff_Ti;dt/coeff_oil,1-dt/coeff_oil];
while err>tol
    i=i+1;
    t(i)=t(i-1)+dt;
    T(:,i)=A*T(:,i-1); %Forward Euler
    err(i)=max([abs(T(1,i)-T(1,i-1))/abs(T(1,i)) abs(T(2,i)-T(2,i-1))/abs(T(2,i))]);
end

%% Grafici

figure(1)
hold on
plot(t,T(2,:)-273.15,'Linewidth',2)
plot(t,T(1,:)-273.15,'Linewidth',2)
xlabel('Tempo [s]')
ylabel('Temperatura [°C]')
legend('Olio','Titanio')
grid minor
set(gca,'Fontsize',16)
print -depsc -f1 ES02Fig1

xx=[0 2500];
figure(2)
hold on
plot(t(2:end),err(2:end),'Linewidth',2)
plot(xx,tol*ones(size(xx)),'g--','Linewidth',2)
xlabel('Tempo [s]')
ylabel('Errore relativo [-]')
grid minor
set(gca,'Fontsize',16)
print -depsc -f2 ES02Fig2

% Calcolo del calore scambiato
q=h*S*(T(1,:)-T(2,:));

figure(3)
hold on
plot(t,q,'Linewidth',2)
xlabel('Tempo [s]')
ylabel('Calore scambiato [W]')
grid minor
set(gca,'Fontsize',16)
print -depsc -f3 ES02Fig3

%% Soluzione numerica (calore specifico funzione della temperatura)
T_2=[T(1,1);T(2,1)];

% Definizione delle funzioni per i calori specifici
cp_oil_fun=@(T)((7200+((T-273.15)/273).^(1.2)));
cp_Ti_fun=@(T)((-468.77*T./(422.11+T).^0.57 - 81.24*T.^2./(1235598+T).^2 + 827.20*T.^3./(51.27+T).^2.79 + 29.98*T.^4./(0.44+T).^3.26));

err=1;
i=1;

dt=10;
t_2=0;

while err>tol
    i=i+1;
    t_2(i)=t_2(i-1)+dt;
    coeff_oil=V_oil*rho_oil*cp_oil_fun(T(2,i-1))/h/S;
    coeff_Ti=V_Ti*rho_Ti*cp_Ti_fun(T(1,i-1))/h/S;
    A=[1-dt/coeff_Ti,dt/coeff_Ti;dt/coeff_oil,1-dt/coeff_oil];
    T_2(:,i)=A*T_2(:,i-1); %Forward Euler
    err=max([abs(T_2(1,i)-T_2(1,i-1))/abs(T_2(1,i)) abs(T_2(2,i)-T_2(2,i-1))/abs(T_2(2,i))]);
end

%% Grafici
% La figura 3 mostra il confronto tra le soluzioni ottenute precedentemente
% (a proprietà costanti) e quelle ottenute con il calore specifico funzione
% della temperatura.
figure(4)
hold on
plot(t,T(2,:)-273.15,'Linewidth',2)
plot(t,T(1,:)-273.15,'Linewidth',2)
plot(t_2,T_2(2,:)-273.15,'--','Linewidth',2)
plot(t_2,T_2(1,:)-273.15,'--','Linewidth',2)
xlabel('Tempo [s]')
ylabel('Temperatura [°C]')
legend('Olio','Titanio','Olio (c_{p,olio} = c_{p,olio}(T))','Titanio (c_{p,Ti} = c_{p,Ti}(T))')
grid minor
set(gca,'Fontsize',16)
print -depsc -f4 ES02Fig4