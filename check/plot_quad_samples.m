close all; clear; clc;


data = load('../data/temp.txt');
figure()
hold on
for i=1:length(data)
plot3(data(i, 1), data(i, 2), data(i, 3), '.k', 'MarkerSize', 20)
end
hold off
axis equal
grid on
grid minor
view([45 45])
