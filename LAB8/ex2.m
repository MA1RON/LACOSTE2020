clc; clear; close all;
%% dati
% geometrici
LL = 1;
DD = 4e-3;
surf = pi*DD^2/4;

% fisici
II = 50; JJ = II/surf;
roel = 1.75e-8;
alpha = 11e-5;
kk = 350;
T0 = 300;
bc_left = load('Tleft.dat');
Tleft = bc_left(:,2);
Tright = T0;

qqq = roel*JJ^2;
rocp = kk/alpha;

%% griglia t
time = bc_left(:,1);
dt = time(2)-time(1);
mm = length(time);

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% BE
Tmid = ones(mm,1)*T0;
aa = dt*alpha/dx^2;
TT = T0*ones(nn,1);
Tp = TT;

diag_si = -aa*ones(nn,1);
diag_p = (1+2*aa)*ones(nn,1);
AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);

% --- condizioni al contorno ---
% dirichlet
AA(1,1) = 1;
AA(1,2) = 0;

% dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;

for inst = 2:mm
    bb = Tp - qqq/kk/rocp*dt;
    % --- condizioni al contorno ---
    % dirichlet
    bb(1) = Tleft(inst);

    % dirichlet
    bb(end) = Tright;
    
    TT = AA\bb;
    Tp = TT;
    
    %{
    plot(xx,TT)
    xlim([0 LL])
    ylim([300 400])
    drawnow
    %}
    
    Tmid(inst) = TT(round(nn/2),1);
end
figure
plot(time/3600,Tmid,'b-','linewidth',2)
grid on
box on
set(gca,'fontsize',18)
xlabel('Tempo [h]')
ylabel('Temperatura centrale [°C]')
title('Transitorio 1D')
