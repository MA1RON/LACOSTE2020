clear; close all; clc;

AA = 1:3;
BB = [4:6]';

try
    CC=[AA;BB];
catch
    CC=[AA;BB'];
end