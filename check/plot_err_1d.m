close all; clear; clc;

x = load('../data/data_x.dat');
y = load('../data/data_y.dat');
data = load('../data/data.dat');

Ns = length(x);

[x, y] = meshgrid(x, y);
x = x';
y = y';

figure()
pcolor(x, y, log10(data))
hold on
plot3([0 1], [0 +1], '-k', 'LineWidth', 2)
plot3([0 1], [0 -1], '-k', 'LineWidth', 2)
hold off
colormap jet
shading flat
axis equal
colorbar
caxis([-6 0])


