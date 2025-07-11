close all; clear; clc;

n = 3;

data_J = load('Fields_J.dat_save');
data_M = load('Fields_M.dat_save');


if n==1

  data = load('data_EJ.txt');

  figure()
  hold on
  plot(data(:, 1), data(:, 2), '-')
  plot(data_J(:, 1), data_J(:, 2), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 3), '-')
  plot(data_J(:, 1), data_J(:, 3), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 4), '-')
  plot(data_J(:, 1), data_J(:, 4), 'o', 'MarkerSize', 4)
  hold off

end

if n==2

  data = load('data_EM.txt');

  figure()
  hold on
  plot(data(:, 1), data(:, 2), '-')
  plot(data_M(:, 1), data_M(:, 2), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 3), '-')
  plot(data_M(:, 1), data_M(:, 3), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 4), '-')
  plot(data_M(:, 1), data_M(:, 4), 'o', 'MarkerSize', 4)
  hold off

end

if n==3
  data = load('data_HJ.txt');

  figure()
  hold on
  plot(data(:, 1), data(:, 2), '-')
  plot(data_J(:, 1), data_J(:, 5), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 3), '-')
  plot(data_J(:, 1), data_J(:, 6), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 4), '-')
  plot(data_J(:, 1), data_J(:, 7), 'o', 'MarkerSize', 4)
  hold off

end

if n==4

  data = load('data_HM.txt');

  figure()
  hold on
  plot(data(:, 1), data(:, 2), '-')
  plot(data_M(:, 1), data_M(:, 5), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 3), '-')
  plot(data_M(:, 1), data_M(:, 6), 'o', 'MarkerSize', 4)
  hold off

  figure()
  hold on
  plot(data(:, 1), data(:, 4), '-')
  plot(data_M(:, 1), data_M(:, 7), 'o', 'MarkerSize', 4)
  hold off

end










