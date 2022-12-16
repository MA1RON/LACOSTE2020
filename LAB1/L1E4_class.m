clear; close all; clc;

MM1 = eye(4);
MM2 = diag(1, 3) + diag(ones(1,3)*10, 1) + diag(logspace(3,9,3),-1) + diag([2 3], -2);
MM3 = MM2';
MM4 = MM1 +1;
MM = [MM1 MM2;MM3 MM4];

%% 4.1
scp = sum(MM(:,2:2:end),2);

%% 4.2
srp = sum(MM(1:2:end,:),1);

%% 4.3
trc = sum(diag(MM));

%% 4.4
trcp = 0;
for j=-6:2:6
    trcp = trcp + sum(diag(MM,j));
end

%% 4.5
MM(MM==1) = -1;

%% 4.6
MM(MM>=0 & MM<=10) = 3;