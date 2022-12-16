clc; clear; close all;
%% dati
% geometrici
DD = .81e-3;
LL = 1;

surf = pi*DD^2/4;

% fisici
alpha = 11e-5;
kk = 350;
rocp = kk/alpha;
T0 = 4.5;
II = 50; JJ = II/surf;
IC = @(T) 50-T.^2;
E0 = 1e-5;
EE = @(T) E0*(II./IC(T)).^1.2;

qqq = @(T) EE(T)*JJ;

%% griglia t
tend = 15e3;
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% BE
Qmid = zeros(mm,1);
Tmid = ones(mm,1)*T0;
TT = T0*ones(nn,1); Tp = TT;
aa = alpha*dt/dx^2;
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
    bb = Tp + dt*qqq(Tp)/rocp; % frozen coefficients
    
    % --- condizioni al contorno ---
    % dirichlet
    bb(1) = T0;

    % dirichlet
    bb(end) = T0;
    
    TT = AA\bb;
    Tp = TT;
    
    %{
    plot(xx,TT)
    ylim([4.5 5.5])
    drawnow
    %}
    
    Tmid(inst) = TT(round(nn/2));
    Qmid(inst) = -kk*(TT(round(nn/2)+1)-2*TT(round(nn/2))+TT(round(nn/2)-1))/dx^2;
end

yyaxis left
plot(time/3600,Tmid,'b-','linewidth',2)
ylabel('Temperatura [K]')
yyaxis right
plot(time/3600,Qmid/1000,'r-','linewidth',2)
ylabel('Potenza [kW/m^3]')
xlabel('Tempo [h]')
grid minor
box on
set(gca,'fontsize',18)
title('Transitorio 1D')