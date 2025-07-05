close all; clear; clc;

x = load('../data/x_bessel.dat');
y = load('../data/y_bessel.dat');
data = load('../data/data_bessel.dat');

Ns = length(x);

[x, y] = meshgrid(x, y);
x = x';
y = y';

figure()
pcolor(x, y, data)
colormap jet
shading flat
axis equal


figure()
pcolor(x, y, abs(besselh(0, 2, x+1j*y)))
colormap jet
shading flat
axis equal
