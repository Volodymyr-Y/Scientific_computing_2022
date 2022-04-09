%% 3D animations
addpath data
addpath figures

%load('test.mat')
load('yy.mat')

% data = soldt01Theat05;
data = soldt1Theat05ImplicitEuler;

rawPressure    = data(:,1:2:end);
rawTemperature = data(:,2:2:end);

numtimesteps   = size(data, 1);
dt             = 1;
t0             = 0;

Nx  = 45;
Ny  = 96;

xx  = linspace(0, 300, Nx);

[gridX, gridY] = meshgrid(xx, yy);

obj = VideoWriter('temp_transient_11.avi');
obj.Quality = 100;
obj.FrameRate = 1/dt;
open(obj);

fig = figure();

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

for i = 1:numtimesteps
    gridP = reshape(rawTemperature(i,:), [Nx, Ny]);
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    
    % fig = figure('position',[10 10 1000 400]);
    
    surface = surf(gridX, gridY, gridP.');
    
    surface.EdgeColor = 'white';
    surface.EdgeAlpha = 0.25;
    
    % view(2)
    view([120 120 120])
    
    xlabel('x-dimension',     'Interpreter','latex')
    ylabel('y-dimension',     'Interpreter','latex')
    zlabel('temperature',     'Interpreter','latex')
    
    title(sprintf('t = %.2f s', t0), 'Interpreter','latex')
    t0 = t0 + dt;

    colormap viridis
    zlim([0 0.8])
    caxis([0 0.8])
    % shading interp;
    colorbar;
    hold off
    
    fname = ['01' num2str(i)]; % full name of image
    print('-djpeg','-r200',fname)     % save image with '-r200' resolution
    I = imread([fname '.jpg']);       % read saved image
    frame = im2frame(I);              % convert image to frame

    writeVideo(obj,frame);

    %     pause(1)
end

obj.close();





%% 2D plots

addpath data
addpath figures

%load('test.mat')
load('yy.mat')

Nx  = 45;
Ny  = 96;

xx  = linspace(0, 300, Nx);

[gridX, gridY] = meshgrid(xx, yy);

fig = figure();

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

gridP = reshape(0*solsteady10(2,:), [Nx, Ny]);
    
set(groot,'defaultAxesTickLabelInterpreter','latex');  

% fig = figure('position',[10 10 1000 400]);

surface = surf(gridX, gridY, gridP.');

surface.EdgeColor = 'white';
surface.EdgeAlpha = 0.25;

view(2)

xlabel('x-dimension',     'Interpreter','latex')
ylabel('y-dimension',     'Interpreter','latex')
zlabel('pressure',     'Interpreter','latex')

colormap viridis

h = colorbar;
set( h, 'YDir', 'reverse' );
set( h, 'TickLabelInterpreter', 'latex');
h.Label.String = 'temperature';
h.Label.Interpreter = 'latex';
h.Label.FontSize = 11;


set(gcf, 'Renderer', 'painters')
print(fig, 'steady_temp_T_0' ,'-dpdf','-r0')

%% Line plot

addpath data
addpath figures

%load('test.mat')
load('yy.mat')

fig = figure();

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

set(groot,'defaultAxesTickLabelInterpreter','latex');  

colors = ["#fde725", "#35b779", "#31688e", "#440154"];

plot(yy(30:end), BL(30:end,1), 'LineWidth', 1.5, 'Color', colors(1))
hold on;
plot(yy(30:end), BL(30:end,2), 'LineWidth', 1.5, 'Color', colors(2))
hold on;
plot(yy(30:end), BL(30:end,3), 'LineWidth', 1.5, 'Color', colors(3))
hold on;
plot(yy(30:end), BL(30:end,4), 'LineWidth', 1.5, 'Color', colors(4))
hold on;
plot(yy(30:end), BL(30:end,4), 'LineWidth', 1.5, 'Color', colors(1), 'LineStyle','--')


xlim([82, 152])

xlabel('y-direction', 'Interpreter', 'latex')
ylabel('temperature', 'Interpreter', 'latex')

legend('$\lambda$ = 0.001', '$\lambda$ = 0.1', '$\lambda$ = 1', '$\lambda$ = 10', '$\lambda$ with Vovas', ...
    'Location', 'best',...
    'Interpreter', 'latex')

% grid on
% 
% set(gcf, 'Renderer', 'painters')
% print(fig, 'boundary_layer_effect' ,'-dpdf','-r0')


%% Mesh plot

load('yy.mat')

Nx  = 45;
Ny  = 96;
zero = zeros(1, Ny*Nx);


xx  = linspace(0, 300, Nx);

[gridX, gridY] = meshgrid(xx, yy);

T = delaunay(gridX, gridY);

TO = triangulation(T,gridX(:), gridY(:), zero(:));

trimesh(TO, 'EdgeColor', '#8E8E8E')
view(2)

ylabel('y-direction', 'Interpreter','latex')
xlabel('x-direction', 'Interpreter','latex')



% 
% set(gcf, 'Renderer', 'painters')
% print(fig, 'mesh_part_1' ,'-dpdf','-r0')








