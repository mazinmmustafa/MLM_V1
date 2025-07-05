close all; clear;

data = load('data.txt');

figure()
hold on
plot(data(:, 1), data(:, 2))
plot(data(:, 1), data(:, 3))
hold off
ylim([0 1])

figure()
hold on
plot(data(:, 1), data(:, 4))
plot(data(:, 1), data(:, 5))
hold off
ylim([0 1])
