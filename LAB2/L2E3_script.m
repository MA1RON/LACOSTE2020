clear; close all; clc;

aa = linspace(-1,1,6);
f1 = L2E3(aa, 1, 2);
f2 = L2E3(aa, 2, 1);

plot(aa, f2, 'linewidth', 3)
hold on 
plot(aa, f1, 'linewidth', 3)
xlabel('a [-]')
ylabel('f [-]')
grid on
box on
set(gca, 'fontsize', 16)

figure()
subplot(1,2,1)
plot(aa, f2, 'linewidth', 3)
subplot(1,2,2)
plot(aa, f1, 'linewidth', 3)