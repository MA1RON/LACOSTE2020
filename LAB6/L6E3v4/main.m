clc; clear; close all;
%% dati
LL = 1;
Text = 300;
hh = 50;

qqqfun = @(x) 8e3*sin(pi*x/LL);
kk = @(T) 100./(11.8+.0238*T) + 8.775e-11*T.^3;

tol = 1e-14;

time = zeros(3,1);
maxjj = 1e5;

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

qqq = qqqfun(xx);

%% NEWTON - derivate numeriche
Tguess = Text*ones(nn,1);

dT = tol;
kder = @(T) (kk(T+dT)-kk(T-dT))/dT;
kder2 = @(T) (kder(T+dT)-kder(T-dT))/dT;

fun = @(T) funbuilder(T,kk,kder,hh,Text,qqq,dx);
jac = @(T) jacbuilder(T,kk,kder,kder2,hh,qqq,dx);

tic
[TT,res,err,lastjj] = newton(fun,jac,Tguess,tol,maxjj);
time(1) = toc;

plot(xx*100,TT+273.15,'r--','linewidth',2,'DisplayName','Newton')
hold on
fprintf('NEWTON numeriche:res = %.2e\terr = %.2e\titerate = %i\n\n',norm(res),err,lastjj)

%% fsolve
Tguess = Text*ones(nn,1);

dT = tol;
kder = @(T) (kk(T+dT)-kk(T-dT))/dT;

fun = @(T) funbuilder(T,kk,kder,hh,Text,qqq,dx);

% OPTIONS = optimoptions('Display','off','TolX',tol,'TolFun',tol);
optimoptions.Display = 'off';

tic
[TT,res,~,OUTPUT] = fsolve(fun,Tguess,optimoptions);
time(2) = toc;

plot(xx*100,TT+273.15,'b-.','linewidth',1,'DisplayName','fsolve')
hold on
fprintf('fsolve:\t\t res = %.2e\titerate = %i\n\n',norm(res),OUTPUT.iterations)

%% punto fisso
% non più fun(TT) = 0 ma TT = fun(TT): fun diversa
Tguess = Text*ones(nn,1);

dT = tol;
kder = @(T) (kk(T+dT)-kk(T-dT))/dT;
kder2 = @(T) (kder(T+dT)-kder(T-dT))/dT;

fun = @(T) funbuilder_fp(T,kk,kder,hh,Text,qqq,dx);

tic
[TT,res,err,lastjj] = punto_fisso(fun,Tguess,tol,maxjj);
time(3) = toc;

plot(xx*100,TT+273.15,'g-.','linewidth',1,'DisplayName','Punto fisso')
fprintf('PUNTO FISSO:\t res = %.2e\terr = %.2e\titerate = %i\n\n',norm(res),err,lastjj)

%% pp
grid on
box on
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('T [^oC]')
legend

figure
bar(categorical({'Newton','fsolve','Punto fisso'}),time)