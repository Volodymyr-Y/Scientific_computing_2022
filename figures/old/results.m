clc; close all; 

Nx = 70;
Ny = 35;

x = linspace(0, 300, Nx);
y = linspace(0, 150, Ny);

[gridX, gridY] = meshgrid(x, y);

gridPvova = reshape(regulargrid7035steadyPressure05, [Nx, Ny]);
gridPvoro = reshape(trigrid7035steadyPressure05, [Nx, Ny]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  

fig = figure('position',[10 10 1000 400]);

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

subplot(1,2,1)
pbaspect([2 1 2])

surface = surf(gridX, gridY, gridPvova.');

surface.EdgeColor = 'white';
surface.EdgeAlpha = 0;

% view(2)
view([120 120 120])

xlabel('x-dimension',     'Interpreter','latex')
ylabel('y-dimension',     'Interpreter','latex')
zlabel('pressure',        'Interpreter','latex')

colormap redblue

shading interp;
colorbar;

pbaspect([2 1 2])

subplot(1,2,2)

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

surface = surf(gridX, gridY, gridPvoro.');

surface.EdgeColor = 'white';
surface.EdgeAlpha = 0;

% view(2)
view([120 120 120])

xlabel('x-dimension',     'Interpreter','latex')
ylabel('y-dimension',     'Interpreter','latex')
zlabel('pressure',        'Interpreter','latex')

colormap redblue

shading interp;
colorbar;

pbaspect([2 1 2])


set(gcf, 'Renderer', 'painters')
print(gcf, '-dpdf', 'steadyPressure0_5.pdf')




