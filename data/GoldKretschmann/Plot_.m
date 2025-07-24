close all; clear; clc;
%%
% data    =   read_gnu_plot_file('data_VED_.txt');
data    =   read_gnu_plot_file('data_HED_.txt');
%%
close all; clc;
figure()
pcolor(data.x, data.y, data.data)
shading flat
colormap jet
colorbar
axis equal
xlabel('$x$ [$\mu$m]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$z$ [$\mu$m]', 'Interpreter', 'Latex', 'FontSize', 14)
title('$20\log_{10}|\mathcal{R}e\{E\}|$ [dBV/m]', 'Interpreter', 'Latex', 'FontSize', 14)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)
set(colorbar, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)
axis([min(min(data.x)) max(max(data.x)) min(min(data.y)) max(max(data.y))])
clim([220 320])
hold on
plot([min(min(data.x)) max(max(data.x))], [0 0], '-k', 'LineWidth', 1)
plot([min(min(data.x)) max(max(data.x))], [-0.05 -0.05], '-k', 'LineWidth', 1)
plot(0, -0.07, 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
hold off

%%
exportgraphics(gcf, 'Figure_.png', 'Resolution', 200)

%%
function[data]=read_gnu_plot_file(filename)
file    =   fopen(filename, 'r');
txt     =   fscanf(file, "%f", [3, inf]);
[~, L]  =   size(txt);
ref     =   txt(2, 1);
for i=2:L
    if ref==txt(2, i)
        break;
    end
end
Ny      =   i-1;
ref     =   txt(1, 1);
for i=2:L
    if ref~=txt(1, i)
        break;
    end
end
Nx      =   i-1;
data.x      =   zeros(Nx, Ny);
data.y      =   zeros(Nx, Ny);
data.data   =   zeros(Nx, Ny);
for i=1:Nx
    for j=1:Ny
        data.x(i, j)    =   txt(1, (i-1)*Nx+i);
        data.y(i, j)    =   txt(2, j);
        data.data(i, j) =   txt(3, (i-1)*Nx+j);
    end
end
fclose(file);
end