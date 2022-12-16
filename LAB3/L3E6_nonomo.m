clc; clear; close all;

%% assegno i dati
tb = .25; %m
ti = .1; %m
kb = .55; %W/mK
ki = .04; %W/mK
Tin = 20; %°C
Text = -5; %°C
hext = 25; %W/m^2K

%% creazione griglie
dx = 5e-5;
xxi = (0:dx:ti)';
nni = length(xxi);

xxb = (ti:dx:tb)'; % sovrappongo un punto in ti!
nnb = length(xxb);

%% assembraggio matrice AA
% ne faccio due diverse con diag diverse per nn diversi
diag_i_i = ones(nni,1);
diag_p_i = -2*ones(nni,1);
diap_s_i = ones(nni,1);
BB_i = [diag_i_i diag_p_i diap_s_i];

AA_i = spdiags(BB_i, -1:1, nni, nni);

diag_i_b = ones(nnb,1);
diag_p_b = -2*ones(nnb,1);
diap_s_b = ones(nnb,1);
BB_b = [diag_i_b diag_p_b diap_s_b];

AA_b = spdiags(BB_b, -1:1, nnb, nnb);


%% assembraggio vettori bb
bb_i = zeros(nni, 1);
bb_b = zeros(nnb, 1);

%% coco
% robin
AA_i(1,1) = -1+hext*dx/-ki; 
AA_i(1,2) = 1; 
bb_i(1) = hext*Text*dx/-ki;
% dirichlet continuità b
AA_b(1,2) = kb; 
AA_b(1,1) = -kb; 
bb_b(1) = ki*(AA_i(end,end)-AA_i(end,end-1));
% dirichlet continuità i
AA_i(end,end-1) = -ki;
AA_i(end,end) = ki; 
bb_i(end) = kb*(AA_b(1,2)-AA_b(1,1));
% dirichlet
AA_b(end,end-1) = 0; 
AA_b(end, end) = 1; 
bb_b(end) = Tin;

%% soluzione numerica
T1 = AA_i\bb_i;
T2 = AA_b\bb_b;
TT = [T1; T2];

%% plotto
plot([xxi; xxb], TT, 'linewidth', 2)
legend('soluzione numerica')
grid on
xlabel('x [m]')
ylabel('T [K]')
