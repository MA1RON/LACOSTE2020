function my_e = L2E4(xx)

my_e = 1;
kk = 1;
while abs((my_e-exp(xx))/exp(xx)) > 1e-14
    my_e = my_e + xx.^kk/factorial(kk);
    kk = kk+1;
end
end