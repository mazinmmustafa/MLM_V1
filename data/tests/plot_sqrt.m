close all; clear; clc;

x = load('data_x.txt');
y = load('data_y.txt');

data_r = load('data_sqrt_r.txt');
data_i = load('data_sqrt_i.txt');


figure()
pcolor(x, y, data_r)
shading flat
colorbar
colormap jet
title('Real')
axis equal

figure()
pcolor(x, y, data_i)
shading flat
colorbar
colormap jet
title('Imagine')
axis equal



