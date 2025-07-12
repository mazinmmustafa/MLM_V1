close all; clear; clc;

data = load('data_PW.txt');
data_E_FEKO = load('Fields_E_PW.dat_save');
data_H_FEKO = load('Fields_H_PW.dat_save');

figure()
hold on
plot(data(:, 1), data(:, 2))
plot(data(:, 1), data(:, 3))
plot(data(:, 1), data(:, 4))
plot(data_E_FEKO(:, 1), data_E_FEKO(:, 2), '--k')
plot(data_E_FEKO(:, 1), data_E_FEKO(:, 3), '--k')
plot(data_E_FEKO(:, 1), data_E_FEKO(:, 4), '--k')
hold off
ylim([0 1.4])

figure()
hold on
plot(data(:, 1), data(:, 5))
plot(data(:, 1), data(:, 6))
plot(data(:, 1), data(:, 7))
plot(data_H_FEKO(:, 1), data_H_FEKO(:, 2), '--k')
plot(data_H_FEKO(:, 1), data_H_FEKO(:, 3), '--k')
plot(data_H_FEKO(:, 1), data_H_FEKO(:, 4), '--k')
hold off
ylim([0 3.5]*1E-3)
