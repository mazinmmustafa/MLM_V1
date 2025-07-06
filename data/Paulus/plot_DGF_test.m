close all; clear; clc;

data = load('data_EJ.txt');

figure()
hold on
plot(data(:, 1)/1E-9, log10(data(:, 2)), '-')
plot(data(:, 1)/1E-9, log10(data(:, 3)), '--')
hold off

figure()
hold on
plot(data(:, 1)/1E-9, log10(data(:, 4)), '-')
plot(data(:, 1)/1E-9, log10(data(:, 5)), '--')
hold off






