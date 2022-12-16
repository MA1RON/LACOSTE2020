clc; clear; close all;
%% dati
% geometrici
LL = 2;
DDint = 0.01;
tt = 1e-3;

surfint = pi*DDint^2/4; surfwa = surfint;
surfext = pi*((DDint/2+tt)^2-DDint^2/4); surfac = surfext;

Awet = pi*DDint*LL;

volint = surfint*LL; volwa = volint;
volext = surfext*LL; volac = volext;
voltot = volint+volext;

% fisici
q0 = 2e3;
qqq = @(x) (q0*exp(-5*(x-1).^2));
GG = 0.12;
rowa = 997;
cpwa = 4186;
roac = 7800;
cpac = 500;
kkac = 30;
hh = 2000;
T0 = 25;

uu = GG/surfint/rowa;

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% griglia t
tend = 20;
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

timetoshow = [dt 1 5 10];
timeshowed = 0;

%% BE - upwind
% solofluido
%{
aa = uu*dt/dx;
qvetwa = qqq(xx)/surfint;
Tp = T0*ones(nn,1);

diag_iwa = (-aa)*ones(size(xx));
diag_pwa = (1+aa)*ones(size(xx));

AA = spdiags([diag_iwa,diag_pwa],[-1 0],nn,nn);
bb = dt*qvetwa/rowa/cpwa;

% --- BC ---
AA(1,1) = 1;
bb(1) = 0;

tol = 1e-10;
res = 1;
inst = 0;

while res > tol
    inst = inst+1;
    TT = AA\(Tp+bb);
    res = abs(TT(end)-Tp(end));
    Tp = TT;
    
    % plot(xx,TT)
    % drawnow
end
%}

% solo tubo
%{
aawa = uu*dt/dx;
betawa = hh*Awet/volint/(rowa*cpwa);
qvetwa = qqq(xx)/surfint;

alpha = kkac/roac/coac;
aaac = alpha*dt/dx^2;
betaac = hh*Awet/volext/(roac*coac);
qqqvetac = qqq(xx)/surfext;

diag_iwa = (-aawa)*ones(size(xx));
diag_pwa = (1+aawa)*ones(size(xx));
diag_iac = (-aaac)*ones(size(xx));
diag_pac = (1+2*aaac+betaac*dt)*ones(size(xx));
diag_sac = diag_iac;

AA = spdiags([diag_iac,diag_pac,diag_sac],[-1 0 1],nn,nn);
bb = (qqqvetac/roac/coac+betaac*T0)*dt;

% --- BC ---
AA(1,1) = 1; 
AA(1,2) = -1;

AA(end,end) = 1;
AA(end,end-1) = -1;

Tp = T0*ones(nn,1);
tol = 1e-5;
res = 1;
inst = 0;

while res>tol
    inst = inst+1;
    BB = Tp+bb;
    
    % --- BC ---
    BB(1) = 0;
    BB(end) = 0;
    
    TT = AA\BB;
    res = abs(max(TT)-max(Tp));
    Tp = TT;
    
    % plot(xx,TT)
    % drawnow
end
%}

% --- risistemo il dominio doppio ---
xxfull=ones(2*nn,1);
xxfull(1:2:end,1)=xx;
xxfull(2:2:end,1)=xx;

% --- proprietà utili per le diagonali ---
aawa = uu*dt/dx;
betawa = hh*Awet/(volwa*rowa*cpwa)*dt;

aaac = kkac/roac/cpac*dt/dx^2;
betaac = hh*Awet/(volac*roac*cpac)*dt;
qqqvetac = qqq(xx)/surfac*dt/(roac*cpac);

% --- diagonali primitive ---
% diag_iwa = (-aawa)*ones(size(xx));
% diag_pwa = (1+aawa)*ones(size(xx));
% diag_iac = (-aaac)*ones(size(xx));
% diag_pac = (1+2*aaac+betaac)*ones(size(xx));
% diag_sac = diag_iac;

% --- diagonali finali ---
%
diag_iifull(1:2:2*nn,1) = (-aaac)*ones(size(xx)); % tubo (i-2)
diag_iifull(2:2:2*nn,1) = (-aawa)*ones(size(xx)); % fluido (i-2)
diag_ifull(1:2:2*nn,1) = -betawa*ones(nn,1); % tubo (i-1) - ! -
diag_ifull(2:2:2*nn,1) = zeros(nn,1); % fluido (i-1) - ! -
diag_pfull(1:2:2*nn,1) = (1+2*aaac+betaac)*ones(size(xx)); % tubo (i)
diag_pfull(2:2:2*nn,1) = (1+aawa+betawa)*ones(size(xx)); % fluido (i)
diag_sfull(1:2:2*nn,1) = zeros(nn,1); % tubo (i+1) - ! -
diag_sfull(2:2:2*nn,1) = -betaac*ones(nn,1); % fluido (i+1) - ! -
diag_ssfull(1:2:2*nn,1) = (-aaac)*ones(size(xx)); % tubo (i+2)
diag_ssfull(2:2:2*nn,1) = zeros(nn,1); % fluido (i+2)
%}
%{
diag_iifull(1:2:2*nn,1) = (-aaac)*ones(size(xx));
diag_iifull(2:2:2*nn,1) = (-aawa)*ones(size(xx));
diag_ifull(1:2:2*nn,1) = -betawa*ones(nn,1); %okkio al conto di spdiags
diag_ifull(2:2:2*nn,1) = zeros(nn,1);
diag_pfull(1:2:2*nn,1) = (1+2*aaac+betaac)*ones(size(xx)); 
diag_pfull(2:2:2*nn,1) = (1+aawa)*ones(size(xx))+betawa;
diag_sfull(1:2:2*nn,1) = zeros(nn,1);
diag_sfull(2:2:2*nn,1) = -betaac*ones(nn,1); %okkio al conto di spdiags
diag_ssfull(1:2:2*nn,1) = (-aaac)*ones(size(xx));
diag_ssfull(2:2:2*nn,1) = zeros(nn,1);
%}
AA = spdiags([diag_iifull, diag_ifull, diag_pfull, diag_sfull ,diag_ssfull],-2:2,2*nn,2*nn);
bb(1:2:2*nn,1) = qqqvetac; % tubo (tn)
bb(2:2:2*nn,1) = zeros(nn,1); % fluido (tn)

% --- condizioni al contorno ---
%{
% tubo iniziale - neumann omogeneo
AA(1,1) = 1;
AA(1,3) = -1;
% bball(1) = 0; % va nel while
% AA(1,2) = 0; % non deve far casino ?

% tubo finale - neumann omogeneo
AA(end-1,end-1) = 1;
AA(end-1,end-3) = -1;
% bball(end) = 0; % va nel while
% AA(end,end) = 1 + aaac; % ???
% AA(end,end-2) = -aaac; % ???
% AA(end-1,end) = 0; % non deve far casino ?
% AA(end,end-1) = 0; % non deve far casino ?

% fluido - dirichlet
AA(2,2) = 1;
AA(2,4) = 0;
% bball(2,2) = T0; % va nel while
% AA(2,1) = 0; % non deve far casino ?
% AA(2,3) = 0; % non deve far casino ?
%}

% tubo iniziale - neumann omogeneo
AA(1,1) = 1;
AA(1,3) = -1;
AA(1,2) = 0;

% fluido iniziale - dirichlet
AA(2,1) = 0;
AA(2,2) = 1;

% tubo finale - neumann omogeneo
AA(end-1,end-1) = 1;
AA(end-1,end) = 0;
AA(end-1,end-3) = -1;

% --- ciclo ---
TTtoshowwa = zeros(nn,length(timetoshow));
TTtoshowac = zeros(nn,length(timetoshow));
TT = T0*ones(2*nn,1);
Tp = T0*ones(2*nn,1);
tol = 1e-4;
res = 1;
inst = 1;
iplot = 1;
while res > tol && inst < mm
    inst = inst+1;
    BB = Tp+bb;
    
    % --- condizioni al contorno ---
    % tubo iniziale - neumann omogeneo
    BB(1) = 0;
    
    % tubo finale - neumann omogeneo
    BB(end-1) = 0;
    
    % fluido iniziale - dirichlet
    BB(2) = T0;
    
    TT = AA\BB;
    TTac = TT(1:2:end,1);
    TTwa = TT(2:2:end,1);
    res = max(abs(max(TTac)-max(Tp(1:2:end))),(TTwa(end)-Tp(end)));
    Tp = TT;
    
    
    figure(1)
    plot(xx,[TTac TTwa],'-','Color',[inst/mm 0 1-inst/mm],'linewidth',2)
    grid on
    box on
    xlabel('x (m)')
    ylabel('Temperature (^oC)')
    set(gca,'fontsize',24)
    title('Transitorio 1D')
    ylim([20 60])
    xlim([0 LL])
    drawnow
    %}
    
    if time(inst) >= timetoshow(timeshowed + 1)
        if timeshowed < length(timetoshow)-1
            timeshowed = timeshowed + 1;
            TTtoshowwa(:,timeshowed) = TTwa;
            TTtoshowac(:,timeshowed) = TTac;
        else
            TTtoshowwa(:,timeshowed + 1) = TTwa;
            TTtoshowac(:,timeshowed + 1) = TTac;
        end
        
    end
end
timetoshow = [timetoshow time(end)];
% salvo steady state
timeshowed = timeshowed + 1;
TTtoshowwa(:,timeshowed + 1) = TTwa;
TTtoshowac(:,timeshowed + 1) = TTac;


figure
for inst = 1:timeshowed + 1
    plot(xx,TTtoshowwa(:,inst),'-','Color',[inst/length(timetoshow) 0 1-inst/length(timetoshow)],'linewidth',2,'DisplayName',['Acqua t=' int2str(timetoshow(inst)) 's'])
    hold on
    plot(xx,TTtoshowac(:,inst),'--','Color',[0 1-inst/length(timetoshow) inst/length(timetoshow)],'linewidth',2,'DisplayName',['Tubo t=' int2str(timetoshow(inst)) 's'])
    hold on
end
legend('location','northeastoutside')
grid on
box on
xlabel('x (m)')
ylabel('Temperature (^oC)')
set(gca,'fontsize',24)
title('Transitorio 1D')
ylim([20 60])