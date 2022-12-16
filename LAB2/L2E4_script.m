clear; close all; clc;

semilogy(readmatrix('valori_esercizio_4.xlsx'), L2E4(readmatrix('valori_esercizio_4.xlsx')), 'ob-.','linewidth',3)
hold on
semilogy(readmatrix('valori_esercizio_4.xlsx'), exp(readmatrix('valori_esercizio_4.xlsx')), 'xr-','linewidth',3)
% semilogy(readmatrix('valori_esercizio_4.xlsx'), (exp(readmatrix('valori_esercizio_4.xlsx'))-L2E4(readmatrix('valori_esercizio_4.xlsx')))/(exp(readmatrix('valori_esercizio_4.xlsx'))), 'dk-.','linewidth',2)