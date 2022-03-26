x = linspace(0, 10, 1000);

figure(1)
hold on;
plot(x,  2*x,  'LineStyle','-',  'Color','k', 'LineWidth', 1.5)
plot(x, -2*x,  'LineStyle','-',  'Color','r', 'LineWidth', 1.5)
plot(x,  x.^2, 'LineStyle','--', 'Color','k', 'LineWidth', 1.5)
plot(x, -x.^2, 'LineStyle','--', 'Color','r', 'LineWidth', 1.5)
hold off;
grid on


figure(2)
x = linspace(0, 300, Nx);
y = linspace(0, 150, Ny);

[gridX, gridY] = meshgrid(x, y);

gridPvova = reshape(regulargrid7035steadyPressure05, [Nx, Ny]);
gridTvova = reshape(regulargrid7035steadyTemperature05, [Nx, Ny]);

gridPvoro = reshape(trigrid7035steadyPressure05, [Nx, Ny]);
gridTvoro = reshape(trigrid7035steadyTemperature05, [Nx, Ny]);

yyaxis right
plot(y, gridPvova(1,:), 'LineStyle','-')
hold on
plot(y, gridPvoro(1,:), 'LineStyle','-')
ylabel('pressure')

yyaxis left
plot(y, gridTvova(1,:), 'LineStyle','-')
hold on
plot(y, gridTvoro(1,:), 'LineStyle','-')
ylabel('pressure')