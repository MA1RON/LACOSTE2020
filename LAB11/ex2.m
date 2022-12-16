clc; clear; close all;
%% dati
% geometrici
LL = 15;
DDint = 15e-3;
DDext = DDint + 2*1.5e-3;
DDiso = DDext + 2*.5e-3;

Awetint = pi*DDint*LL;
Awetext = pi*DDext*LL;

surfint = pi*DDint^2/4;
surfext = pi*(DDext^2-DDint^2)/4;
surfiso = pi*(DDiso^2-DDext^2)/4;

volint = surfint*LL; % water
volext = surfext*LL; % copper
voliso = surfiso*LL; % insulating
voltu = volext + voliso; % tube

% fisici
rowa = 997;
cpwa = 4168;
muwa = 1.12e-6;
kkwa = .6;

kkcu = 385.0;
rocu = 8960;
cpcu = 385;

kkiso = .1;

kktu = (kkcu*volext + kkiso*voliso)/voltu;
% rotube ~ rocu
% cptube ~ cpcu

II = 1e3;
GG = .2; uu = GG/rowa/surfint; 
ReDD = uu*DDint*rowa/muwa;
Pr = cpwa*muwa/kkwa;
pp = 2e5;
TTair = 20;
T0 = 5;

roel = 1.68e-8;
qqq = roel*(II/surfext)^2;

NuDD = .023*ReDD^(4/5)*Pr^(.4); % fonte: Tabelle calcoli termici e idraulici
hhint = NuDD*kkwa/DDint;
hhext = kkiso/DDiso*(.6)^2; % fonte: Wikipedia, Churchill and Chu correlation

%% griglia spaziale
% fluido 1 - 3 - 5 - ...
% tubo + isolante 2 - 4 - 6 - ...
dx = LL/2/100;
xx = (0:dx:LL)';
xx = sort([xx;xx]);
nn = length(xx);

%% griglia temporale
tend = 70;
dt = tend/50;
time = (0:dt:tend)';
mm = length(time);

%% sistema lineare (BE + upwind)
% parametri utili
aawa = dt*uu/dx;
betawa = dt*hhint*Awetint/rowa/volint/cpwa;
aatu = dt*kktu/rocu/cpcu/dx^2;
qqqtu = qqq*dt/rocu/cpcu;
betatuint = dt*hhint*Awetint/voltu/rocu/cpcu;
betatuext = dt*hhext*Awetext/voltu/rocu/cpcu;

% diagonali
diag_ii(1:2:2*nn,1) = -aawa*ones(nn,1); % fluido (i-2)
diag_ii(2:2:2*nn,1) = -aatu*ones(nn,1); % tubo (i-2)
diag_i(1:2:2*nn,1) = -betatuint*ones(nn,1); % fluido (i-1)
diag_i(2:2:2*nn,1) = zeros(nn,1); % tubo (i-1)
diag_p(1:2:2*nn,1) = (1 + aawa + betawa)*ones(nn,1); % fluido (i)
diag_p(2:2:2*nn,1) = (1 + 2*aatu + betatuint + betatuext)*ones(nn,1); % tubo (i)
diag_s(1:2:2*nn,1) = zeros(nn,1); % fluido (i+1)
diag_s(2:2:2*nn,1) = -betawa*ones(nn,1); % tubo (i+1)
diag_ss(1:2:2*nn,1) = zeros(nn,1); % fluido (i+2)
diag_ss(2:2:2*nn,1) = -aatu*ones(nn,1); % tubo (i+2)

% termini matriciali
TT = T0*ones(2*nn,1); Tp = TT;
TTairvet = TTair*ones(nn,1);

AA = spdiags([diag_ii diag_i diag_p diag_s diag_ss],-2:2,2*nn,2*nn);
bb(1:2:2*nn,1) = zeros(nn,1); % fluido (i)
bb(2:2:2*nn,1) = (qqqtu + betatuext*TTair)*ones(nn,1); % tubo (i)

% condizioni al contorno
% fluido iniziale - dirichlet
AA(1,1) = 1;
AA(1,2) = 0;

% fluido finale - come la matrice
% AA(end-1,end) = 0;
% AA(end-1,end-1) = 1;
% AA(end-1,end-2) = 0;
% AA(end-1,end-3) = -1;

% tubo iniziale - neumann omogeneo
AA(2,1) = 0;
AA(2,2) = 1;
AA(2,3) = 0;
AA(2,4) = -1;

% tubo finale - neumann omogeneo
AA(end,end-2) = -1;
AA(end,end-1) = 0;
AA(end,end) = 1;

% ciclo
for inst = 2:mm
    BB = Tp + bb;
    
    % condizioni al contorno
    BB(1) = T0; % fluido iniziale - dirichlet
    % BB(end-1) = 0; % fluido finale - come la matrice
    BB(2) = 0; % tubo iniziale - neumann omogeneo
    BB(end) = 0; % tubo finale - neumann omogeneo
    
    TT = AA\BB;
    
    TTwa = TT(1:2:end);
    TTtu = TT(2:2:end);
    
    plot(xx,[TTwa TTtu TTairvet],'linewidth',3)
    legend('Acqua','Tubo','Aria','location','northeastoutside')
    grid on
    box on
    xlabel('x (m)')
    ylabel('Temperature (^oC)')
    set(gca,'fontsize',24)
    title('Transitorio 1D')
    ylim([0 21])
    xlim([0 LL])
    drawnow
    
    Tp = TT;
end