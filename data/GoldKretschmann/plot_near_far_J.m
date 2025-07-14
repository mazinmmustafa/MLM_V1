close all; clear; clc;

filename = 'data_near_far_J';

data = load(strcat(filename, '.txt'));

E_theta = data(:, 2);
E_theta(isinf(E_theta)) = -1.0E16;

dB_min = +220;

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

E_theta1 = E_theta;

data_save = [(data(:, 1))  E_theta];

figure()
polar(flipud(data(:, 1))+pi/2, E_theta)

E_theta = data(:, 3);

E_theta(isinf(E_theta)) = -1.0E16;


E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

data_save = [data_save E_theta];

E_theta2 = E_theta;

figure()
hold on
polar(flipud(data(:, 1))+pi/2, E_theta1)
polar(flipud(data(:, 1))+pi/2, E_theta2, '--')
hold off
axis equal
##

file = fopen(strcat(filename,'_.txt'), 'w');
Theta = linspace(0, 2*pi, length(E_theta));
for i=1:length(data_save)
  fprintf(file, '%21.14E %21.14E %21.14E \n',...
  data_save(i, 1), data_save(i, 2)+dB_min, data_save(i, 3)+dB_min);
end
fclose(file);



