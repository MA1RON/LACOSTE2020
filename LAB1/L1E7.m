clear; close all; clc;

nn = 10;
FF = zeros(nn,1);
FF(2) = 1;

for jj=3:nn
    FF(jj) = FF(jj-1) + FF(jj-2);
end