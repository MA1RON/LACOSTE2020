clear; close all; clc;

%% metodo 1
tic
for jj = 1:1e6
    v1(jj) = jj;
end
t1 = toc;

%% metodo 2
tic
v2 = zeros(1,1e6);
for jj = 1:1e6
    v2(jj) = jj;
end
t2 = toc;

%% metodo 3
tic
v3 = 1:1e6;
t3 = toc;