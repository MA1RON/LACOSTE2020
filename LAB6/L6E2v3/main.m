clc; clear; close all;
%% dati
ro = 999;
mu = 1.12e-3;
gg = 9.805;
LL = 152;
DD = .23;
KK = .8+5+1.5*4+1;
hl = 76-60;

tol = 1e-12;

% in: in(1) = f; in(2) = v
fun1 = @(f,v) f - (.79*log(ro*v*DD/mu) - 1.64).^-2;
fun2 = @(f,v) (f*LL/DD + KK) .* v.^2 - 2*gg*hl;

fun = @(in) [fun1(in(1),in(2)); fun2(in(1),in(2))];

% derivate numeriche
df = tol; dv = tol;
dfun1dv = @(f,v) (fun1(f,v+dv)-fun1(f,v-dv))/(2*dv);
dfun1df = @(f,v) (fun1(f+df,v)-fun1(f-df,v))/(2*df);
dfun2dv = @(f,v) (fun2(f,v+dv)-fun2(f,v-dv))/(2*dv);
dfun2df = @(f,v) (fun2(f+df,v)-fun2(f-df,v))/(2*df);

jac = @(in) [dfun1df(in(1),in(2)),dfun1dv(in(1),in(2));...
             dfun2df(in(1),in(2)),dfun2dv(in(1),in(2))];
         
%% NEWTON
in0 = [.001; .35];

tic
[out,res,err,lastjj] = newton(fun,jac,in0,tol,1e3);
fprintf('NEWTON NUMERICO: f = %.4f\t\tv = %.4fm/s\n',out(1),out(2))
fprintf('\t\t err = %.3e\tIterate = %i\t\tResiduo = %.3e\n',err,lastjj,norm(res))
fprintf('\t\t tempo = %.3f s\tTolleranza = %.3e\n\n',toc,tol)