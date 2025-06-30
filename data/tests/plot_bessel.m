close all; clear; clc;

x = load('data_x.txt');
y = load('data_y.txt');

data_m = load('data_besselh_m.txt');
data_p = load('data_besselh_p.txt');


figure()
pcolor(x, y, log10(data_m))
shading flat
colorbar
colormap jet

figure()
pcolor(x, y, data_p)
shading flat
colorbar
colormap jet

z = x+1j*y;
ans = besselh(0, 2, z);

figure()
pcolor(x, y, log10(abs(ans)))
shading flat
colorbar
colormap jet

figure()
pcolor(x, y, angle(ans))
shading flat
colorbar
colormap jet


