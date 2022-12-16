clear; close all; clc;

tt = 10:10:50;
c = categorical({['t_1 = ' num2str(tt(1)) ' [s]'],...
    ['t_2 = ' num2str(tt(2)) ' [s]'],...
    ['t_3 = ' num2str(tt(3)) ' [s]'],...
    ['t_4 = ' num2str(tt(4)) ' [s]'],...
    ['t_5 = ' num2str(tt(5)) ' [s]']});
TT = readmatrix('valori_esercizio_7.xlsx');

h=bar(c, TT,.4);
h.FaceColor = 'flat';
h.CData = [0 0 1];
hold on
plot(c, TT,'or-', 'linewidth',3)
ylim([275 305])
ylabel('Temperatura [K]')
set(gca, 'fontsize', 20)
grid on