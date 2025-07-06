close all; clear; clc;

data = load('data_v.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
plot(data(:, 1), data(:, 3), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
plot(data(:, 1), data(:, 5), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 6), '-')
plot(data(:, 1), data(:, 7), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 8), '-')
plot(data(:, 1), data(:, 9), '--')
hold off

data = load('data_i.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
plot(data(:, 1), data(:, 3), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
plot(data(:, 1), data(:, 5), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 6), '-')
plot(data(:, 1), data(:, 7), '--')
hold off

figure()
hold on
plot(data(:, 1), data(:, 8), '-')
plot(data(:, 1), data(:, 9), '--')
hold off
