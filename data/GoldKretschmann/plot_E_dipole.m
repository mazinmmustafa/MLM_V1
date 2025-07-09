close all; clear; clc;

x = load('data_x.txt');
z = load('data_z.txt');

data = load('data.txt');

figure()
pcolor(x/1E-6, z/1E-6, 20*log10(data))
hold on
plot([-2 +2], 0*[1 1], '-k', 'LineWidth', 1)
plot([-2 +2], -0.05*[1 1], '-k', 'LineWidth', 1)
hold off
shading flat
axis equal
colorbar
colormap jet
caxis([220 320])
