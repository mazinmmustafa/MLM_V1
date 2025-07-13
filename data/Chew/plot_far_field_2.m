close all; clear; clc;

filename = 'data_far_field_J_2';

data = load(strcat(filename, '.txt'));

E_theta = data(:, 2);
E_theta(isinf(E_theta)) = -1.0E16;

dB_min = -20;

E_max = max(10.^(E_theta/20));

E_theta = E_theta-20*log10(E_max);

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

E_theta1 = E_theta;

data_save = [(data(:, 1))  E_theta];

figure()
polar(flipud(data(:, 1))+pi/2, E_theta)

E_theta = data(:, 3);

E_theta(isinf(E_theta)) = -1.0E16;


E_max = max(10.^(E_theta/20));
E_theta = E_theta-20*log10(E_max);

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

data_save = [data_save E_theta];

E_theta2 = E_theta;

figure()
hold on
polar(flipud(data(:, 1))+pi/2, E_theta1)
polar(flipud(data(:, 1))+pi/2, E_theta2)
hold off
axis equal



