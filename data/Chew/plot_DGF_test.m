close all; clear; clc;

data = load('data_EJ.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 3), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
hold off

data = load('data_EM.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 3), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
hold off

data = load('data_HJ.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 3), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
hold off

data = load('data_HM.txt');

figure()
hold on
plot(data(:, 1), data(:, 2), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 3), '-')
hold off

figure()
hold on
plot(data(:, 1), data(:, 4), '-')
hold off









