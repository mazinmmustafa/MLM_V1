close all; clear; clc;

data_EJ = load('data_far_field_EJ.txt');

figure()
hold on
plot(data_EJ(:, 1), 20*log10(data_EJ(:, 2)))
plot(data_EJ(:, 1), 20*log10(data_EJ(:, 3)))
hold off
