clear; close all; clc;

QQ = diag([3 4], -2) + diag(pi/2:pi/2:3*pi/2, -1) + diag(-12:5:3) + diag(3*pi/2:pi/2:5*pi/2, 1) + diag([-2 -1], 2);
QQ(end) = 1;

%% 6.1
for ii = 1:size(QQ,1)
    for jj = 1:size(QQ,2)
        if ii==jj
            QQ(ii,jj) = 10;
        end
    end
end

%% 6.2
while QQ(2,1)<QQ(1,1)
    QQ(2,1) = QQ(2,1) + 1;
end

%% 6.3
trcp = 0;
for ii = 1:size(QQ,1)
    for jj = 1:size(QQ,2)
        if ii==jj
            trcp = trcp + QQ(ii,jj);
        end
    end
end

%% 6.4
QQ(~QQ) = 1;
QQ(QQ == 0) = 1;
%{
for ii1 = 1:size(QQ,1)
    for jj1 = 1:size(QQ,2)
        for ii2 = 1:size(QQ,1)
            for jj2 = 1:size(QQ,2)
                if QQ(ii1,jj1) == QQ(ii2, jj2) &&~(ii1 == ii2 && jj1 == jj2)
                    QQ(ii1,jj1) = 0;
                    QQ(ii2, jj2) = 0;
                end
            end
        end
    end
end
%}
%{
for ii = 1:size(QQ,1)
    for jj = 1:size(QQ,2)
        if length(find(QQ == QQ(ii,jj))) > 1 && QQ(ii,jj) ~= 0
            QQ(QQ==QQ(ii,jj)) = 0;
        end
    end
end
%}