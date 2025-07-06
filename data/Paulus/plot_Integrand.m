close all; clear; clc;

data = load('data_I.txt');

figure()
hold on
plot(data(:, 1), abs(data(:, 2)))
hold off
ylim([0 1.5E10])

figure()
hold on
plot(data(:, 1), abs(data(:, 3)))
hold off
ylim([0 1.5E10])

figure()
hold on
plot(data(:, 1), abs(data(:, 4)))
hold off
ylim([0 3.0E8])

