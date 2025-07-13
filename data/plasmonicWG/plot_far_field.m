close all; clear; clc;

filename = 'data_far_field.txt';

data = load(filename);

E_theta = data(:, 2);
E_theta(isinf(E_theta)) = -1.0E10;

dB_min = -40;

E_max = max(10.^(E_theta/20))
E_theta = E_theta-20*log10(E_max);

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

data_save = [(data(:, 1))  E_theta];

figure()
polar(flipud(data(:, 1))+pi/2, E_theta)

E_theta = data(:, 3);

E_theta(isinf(E_theta)) = -1.0E10;

dB_min = -40;

E_max = max(10.^(E_theta/20))
E_theta = E_theta-20*log10(E_max);

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

data_save = [data_save E_theta];

figure()
polar(flipud(data(:, 1))+pi/2, E_theta)


file = fopen(filename, 'w');
Theta = linspace(0, 2*pi, length(E_theta));
for i=1:length(data_save)
  fprintf(file, '%21.14E %21.14E %21.14E\n',...
  data_save(i, 1), data_save(i, 2)+dB_min, data_save(i, 3)+dB_min);
end
fclose(file);


