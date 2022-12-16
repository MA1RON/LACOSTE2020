clear; close all; clc;

nn = 10; kk = 1e5;
AA = rand(nn/2,nn); BB = randi(kk,nn/2,nn)/kk;
CC = zeros(nn,nn); 
CC(1:nn/2,:) = AA; 
CC(nn/2+1:end,:) = BB;

%% 10.1
vv = find(CC > 1e-3);
xx = zeros(length(vv),1);
yy = zeros(length(vv),1);
for jj = 1:length(vv)
    yy(jj) = ceil(vv(jj)/nn);
    xx(jj) = rem(vv(jj),nn);
end


%% 10.2
CC = ones(nn);
cont = 0;
while length(find(CC > 1e-3)) == nn^2
    AA = rand(nn/2,nn); BB = randi(kk,nn/2,nn)/kk; 
    CC(1:nn/2,:) = AA; CC(nn/2+1:end,:) = BB;
    cont = cont+1;
end

%% 10.3
AA = rand(nn/2,nn); BB = randi(kk,nn/2,nn)/kk; 
CC(1:nn/2,:) = AA; CC(nn/2+1:end,:) = BB;
CC(:,2:2:end) = 0;
while sum(CC,'all')>10
    if sum(CC,'all') > 30
        CC(1:2:end,:) = CC(1:2:end,:) -1;
        CC(:,2:2:end) = CC(:,2:2:end) -2;
    elseif sum(CC,'all') > 25
        CC(1,:) = CC(1,:)/2;
    else
        CC = CC/2;
    end
end
