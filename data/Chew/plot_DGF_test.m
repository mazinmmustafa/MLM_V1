close all; clear; clc;

n = 1;

if n==1
  
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

end

if n==2
  
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

end

if n==3
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

end

if n==4
  
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
  
end










