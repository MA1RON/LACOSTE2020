clc; clear; close all;
%% dati
kk = 1.3806503e-23;
aa = .401;
bb = 42.7e-6;
NN = 1e3;
TT = 300;
pp = 3.5e7;

fun = @(v) (pp+aa*(NN./v).^2).*(v-NN*bb) - kk*NN*TT;

tol = 1e-12;

v0 = .01; % guess iniziale

%% NEWTON - derivata analitica
der = @(v) (-2*aa*NN^2./(v.^3)).*(v-NN*bb) + pp + aa*(NN./v).^2;

tic
[vv,res,err,lastjj] = newton(fun,der,v0,tol,1e3);
fprintf('NEWTON derivata analitica: v =  %.5f m^3\n',vv);
fprintf('res = %.2e\terr = %.2e\titerate = %.2i\ttempo = %.4fs\n\n',res,err,lastjj,toc);

%% NEWTON - derivata numerica
dv = tol;
der = @(v) (fun(v+dv)-fun(v-dv))/2/dv;

tic
[vv,res,err,lastjj] = newton(fun,der,v0,tol,1e3);
fprintf('NEWTON derivata numerica: v =  %.5f m^3\n',vv);
fprintf('res = %.2e\terr = %.2e\titerate = %.2i\ttempo = %.4fs\n\n',res,err,lastjj,toc);

%% BISEZIONE
v0_bisez = [v0 v0*100];

tic
[vv,res,err,lastjj] = bisezione(fun,v0_bisez,tol,1e3);
fprintf('BISEZIONE: v =  %.5f m^3\n',vv);
fprintf('res = %.2e\terr = %.2e\titerate = %.2i\ttempo = %.4fs\n\n',res,err,lastjj,toc);

%% fzero
OPTIONS = optimset('TolX',tol,'TolFun',tol);

tic
[vv,res,~,OUTPUT] = fzero(fun,v0,OPTIONS);
fprintf('fzero: v =  %.5f m^3\n',vv);
fprintf('res = %.2e\titerate = %.2i\ttempo = %.4fs\n\n',res,OUTPUT.iterations,toc);