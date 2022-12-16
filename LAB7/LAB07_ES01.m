%% LaCoSTe 2020/2021 
%% LAB07_ES01 - Soluzioni, Matteo Nicoli
clear all
close all
clc

%% Dati di input
%Temperature iniziali
T=450;
Tw=15;

% Scambio termico
h=100;
S=2*pi/4*(0.6^2)+2*pi*0.4*0.09+pi*0.6*0.12;

% Proprietà alluminio
V_al=pi/4*(0.4^2*(0.3-0.1)+0.6^2*0.1);
k_al=237;
rho_al=2700;
c_al=897;

% Proprietà acciaio
V_st=pi/4*((0.6^2-0.4^2)*0.01)*2;
k_st=14.2;
rho_st=7978;
c_st=480;

% Proprietà parte meccanica
V=V_al+V_st;
k=(k_al*V_al+k_st*V_st)/V;
rho=(rho_al*V_al+rho_st*V_st)/V;
c=(c_al*rho_al*V_al+c_st*rho_st*V_st)/(rho*V);

% Numero di Biot e tau
L=V/S;
Bi=h*L/k;
tau=V*rho*c/h/S;

%% Soluzione numerica
tol=1e-4;
err=1;
i=1;

% Imposto il dt e l'istante iniziale
dt=50;
t=0;

coeff=V*rho*c/h/S;
while err(i)>tol
    i=i+1;
    T(i)=T(i-1)+dt/coeff*(Tw-T(i-1)); %Forward Euler
    %T(i)=1/(1+dt/coeff)*(T(i-1)+dt/coeff*Tw); %Backward Euler
    t(i)=t(i-1)+dt;
    err(i)=abs(T(i)-T(i-1))/abs(T(i));
end

%% Soluzione analitica
T_an=@(t)((T(1)-Tw)*exp(-t/coeff)+Tw);

%% Grafici
figure(1)
hold on
plot(t,T,'Linewidth',2)
plot(t,T_an(t),'--','Linewidth',2)
xlabel('Tempo [s]')
ylabel('Temperatura [°C]')
legend('Soluzione numerica','Soluzione analitica')
grid minor
set(gca,'Fontsize',16)
print -depsc -f1 ES01Fig1

xx=0:15000;
figure(2)
hold on
plot(t(2:end),err(2:end),'Linewidth',2)
plot(xx,tol*ones(size(xx)),'g--','Linewidth',2)
xlabel('Tempo [s]')
ylabel('Errore relativo [-]')
grid minor
set(gca,'Fontsize',16)
print -depsc -f2 ES01Fig2