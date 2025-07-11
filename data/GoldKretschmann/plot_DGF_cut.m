close all; clear; clc;

data = load('data_cut.txt');

figure()
hold on
plot(data(:, 1), 20*log10(data(:, 4)))
hold off
##ylim([200 400])