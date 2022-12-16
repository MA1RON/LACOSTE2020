clear; close all; clc;

%% 11.1
nn = 10;
coo = rand(nn,2);

%% 11.2
mm = 0;
for cc=1:nn
    if norm([coo(cc,1) coo(cc,2)]) <= 1
        mm = mm+1;
    end
end

%% 11.3
my_pi = 4*mm/nn;

%% 11.4
my_pi = 0;
nn = 1e4;
err = ones(nn);
while err(log10(nn))>=1e-3 % 1e-6 non ce la fa
    nn = nn*10;
    xx = rand(nn,1);
    yy = rand(nn,1);
    mm = sum(sqrt(xx(:).^2+yy(:).^2) <= 1);

    my_pi = 4*mm/nn;
    err(log10(nn)) = abs(my_pi-pi)/pi;
end

loglog(logspace(1,log10(nn),log10(nn)), err)