
clear; clc; close all; 

load meshIndependency.mat
fn = fieldnames(meshIndependency);
vector      = zeros(numel(fn)-1, 1);
numelvector = zeros(numel(fn)-1, 1);

for i=1: numel(fn)
    
    Nx = meshIndependency.(fn{i}).nx;
    Ny = meshIndependency.(fn{i}).ny;
    
    x = linspace(0, 300, Nx);
    y = linspace(0, 150, Ny);
    
    [gridX, gridY] = meshgrid(x, y);
    
    gridP = reshape(meshIndependency.(fn{i}).pres, [Nx, Ny]);
    
    
    % Plotting 
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    fig = figure(i);
    
    set(fig, 'Units', 'centimeters')
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    
    surface = surf(gridX, gridY, gridP.');
    
    surface.EdgeColor = 'white';
    surface.EdgeAlpha = 0;
    
    colormap turbo
    
    % view(2)
    view([120 120 120])
    
    xlabel('x-dimension',     'Interpreter','latex')
    ylabel('y-dimension',     'Interpreter','latex')
    zlabel('pressure',        'Interpreter','latex')
    % 
    % set(gcf, 'Renderer', 'painters')
    % print(gcf, '-dpdf', 'steadyTemperature2Dplot.pdf')

    meshIndependency.(fn{i}).averageP = mean(meshIndependency.(fn{i}).pres);

    if i > 1
        meshIndependency.(fn{i}).averagePerror = meshIndependency.(fn{i}).averageP - meshIndependency.(fn{i-1}).averageP;
        vector(i-1) = meshIndependency.(fn{i}).averagePerror;
        numelvector(i-1) = Nx * Ny;
    end

end
