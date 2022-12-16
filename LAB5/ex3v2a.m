clc; clear; close all;
%% dati
% dominio
% Tint | 0 - ss - contatto - is - end | Text
% geometrici
tss = 1e-3;
tis = 2e-3;

% fisici
Tint = 70;
Text = 20;
hhext = 15;
kkss = 15;
kkis = 1;

%% griglia
dx = 1e-4;
xx = [(0:dx:tss)';(tss+dx:dx:tss+tis)'];
nncontatto = length(0:dx:tss);
nn = length(xx);

%% costruisco il sistema
diag_si = ones(nn,1);
diag_p = -2*ones(nn,1);
AA = spdiags([diag_si diag_p diag_si],-1:1,nn,nn);
bb = zeros(nn,1);

% --- condizioni al contorno ---
% dirichlet interno
AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Tint;

% contatto
AA(nncontatto,nncontatto-1) = -1;
AA(nncontatto,nncontatto) = 1+kkis/kkss;
AA(nncontatto,nncontatto+1) = -kkis/kkss;
bb(nncontatto) = 0; % pleonastico

% robin esterno
AA(end,end-1) = 1;
AA(end,end) = -1-hhext*dx/kkis;
bb(end) = -hhext*Text*dx/kkis;

%% risolvo il sistema
TT = AA\bb;

%% post production
plot(xx*1e3,TT,'b-.','linewidth',2)
hold on
plot([xx(nncontatto) xx(nncontatto)]*1e3,[min(TT) max(TT)],'k--','linewidth',2)
legend('Temperatura', 'Contatto', 'location', 'northeast')
xlabel('Spessore [mm]')
ylabel('Temperatura [°C]')
grid on
box on
set(gca, 'fontsize', 18)