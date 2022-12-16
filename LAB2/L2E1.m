clear; close all; clc;

xx = linspace(-2,2);
f = @(x)1./(1+x.^2);
plot(xx,f(xx),'r-.','linewidth',3)
xlabel('X [-]')
ylabel('1/(1+X^2) [-]')
legend('f(x)','location','north')
grid on % minor
box on
set(gca, 'fontsize', 16)
print('-f1', '-depsc', 'L2E1_eps')
print('-f1', '-djpeg', 'L2E1_jpeg')