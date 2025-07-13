close all; clear; clc;

filename = 'data_far_field_J';
##filename = 'data_far_field_M';

data = load(strcat(filename, '.txt'));

E_theta = data(:, 2);
E_theta(isinf(E_theta)) = -1.0E16;

dB_min = -20;

E_max = max(10.^(E_theta/20));

E_theta = E_theta-20*log10(E_max);

E_theta(E_theta<dB_min) = dB_min;
E_theta = E_theta-dB_min;

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

figure()
polar(flipud(data(:, 1))+pi/2, E_theta)


E_phi = data(:, 4);
E_phi(isinf(E_phi)) = -1.0E16;

dB_min = -20;

E_max = max(10.^(E_phi/20));
E_phi = E_phi-20*log10(E_max);

E_phi(E_phi<dB_min) = dB_min;
E_phi = E_phi-dB_min;

data_save = [data_save  E_phi];

figure()
polar(flipud(data(:, 1))+pi/2, E_phi)


E_phi = data(:, 5);
E_phi(isinf(E_phi)) = -1.0E16;


E_max = max(10.^(E_phi/20));
E_phi = E_phi-20*log10(E_max);

E_phi(E_phi<dB_min) = dB_min;
E_phi = E_phi-dB_min;

data_save = [data_save  E_phi];

figure()
polar(flipud(data(:, 1))+pi/2, E_phi)


file = fopen(strcat(filename,'_.txt'), 'w');
Theta = linspace(0, 2*pi, length(E_theta));
for i=1:length(data_save)
  fprintf(file, '%21.14E %21.14E %21.14E %21.14E %21.14E\n',...
  data_save(i, 1), data_save(i, 2)+dB_min, data_save(i, 3)+dB_min,...
   data_save(i, 4)+dB_min, data_save(i, 5)+dB_min);
end
fclose(file);


