clc; clear; close all;
%% dati
% geometrici
DD = 12e-3;

% fisici
kk = 40;
ro = 7800; 
cc = 600;
alpha = kk/ro/cc;
Tinit = 1150;
Tfinal = 400;
Tair = @(t) 325 + 0.1875*t;
hh = 20;

%% griglia r
dr = DD/2/100;
rr = (0:dr:DD/2)';
nn = length(rr);

%% griglia t
tend = 10000; % guess
dt = tend/1000;
time = (0:dt:tend)';
mm = length(time);

%% BE
TTave = ones(size(time))*Tinit;
TTsup = TTave;
TT = Tinit*ones(nn,1);
Tp = TT;

aa = alpha*dt;
diag_s = - aa/2./rr/dr - aa/dr^2;
diag_p = (1 + 2*aa/dr^2)*ones(nn,1);
diag_i = + aa/2./rr/dr - aa/dr^2;
AA = spdiags([[diag_i(2:end);0] diag_p [0;diag_s(1:end-1)]],-1:1,nn,nn);

% --- condizioni al contorno ---
% neumann
AA(1,1) = 1;
AA(1,2) = -1;

% robin
AA(end,end-1) = 1;
AA(end,end) = - 1 - hh/kk*dr;

inst = 2;
while inst < mm && TTave(inst-1) > Tair(time(inst-1))
    bb = Tp;
    
    % --- condizioni al contorno ---
    % neumann
    bb(1) = 0;

    % robin
    bb(end) = -hh/kk*dr*Tair(time(inst));
    
    TT = AA\bb;
    
    Tp = TT;
    
    %{
    plot(rr,TT)
    ylim([0 3000])
    drawnow
    %}
    TTave(inst) = sum(TT)/nn;
    TTsup(inst) = TT(end);
    
    inst = inst +1;
end

plot(time(1:inst-1),TTave(1:inst-1),'b-','linewidth',2)
hold on
plot(time(1:inst-1),TTsup(1:inst-1),'r-.','linewidth',2)
hold on
plot(time(inst-1),TTave(inst-1),'ko','markerfacecolor','y','markersize',10)
hold on
plot(time(1:inst-1),Tair(time(1:inst-1)),'k--','linewidth',2)
legend('T_{media}','T_{superficie}','T_{minima}','T_{aria}','location','northeast')
grid on
box on
xlabel('Tempo [s]')
ylabel('Temperatura [K]')
set(gca,'fontsize',18)
title('Transitorio 1D')