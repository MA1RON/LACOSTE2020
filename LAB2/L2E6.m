clear; close all; clc;

xx = linspace(0,10);
f1 = @(x) x.^2 - x + 1;
f2 = @(x) exp(x);

yyaxis left
plot(xx,f1(xx), 'b', 'linewidth', 3)
ylabel('Temperatura [°C]')

yyaxis right
plot(xx,f2(xx), 'r-.', 'linewidth', 3)
ylabel('Potenza Termica [W]')

xlabel('x [mm]')
grid on
legend('Temperatura [°C]', 'Potenza Termica [W]', 'location', 'northwest')
set(gca, 'fontsize', 16)