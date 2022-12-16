function T = L2E5(X)
T = X;
for ii=1:size(X,1)
    for jj=1:size(X,2)
        T(ii,jj) = X(jj,ii);
    end
end
end
