%% 3D animations
addpath data
addpath figures

%load('test.mat')
load('yy.mat')

% data = soldt01Theat05;
data = soldt1Theat10;

rawPressure    = data(:,1:2:end);
rawTemperature = data(:,2:2:end);

numtimesteps   = size(data, 1);
dt             = 1;
t0             = 0;

Nx  = 45;
Ny  = 96;

xx  = linspace(0, 300, Nx);

[gridX, gridY] = meshgrid(xx, yy);

obj = VideoWriter('press_transient_1.avi');
obj.Quality = 100;
obj.FrameRate = 1/dt;
open(obj);

fig = figure();

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

for i = 1:numtimesteps
    gridP = reshape(rawPressure(i,:), [Nx, Ny]);
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    
    % fig = figure('position',[10 10 1000 400]);
    
    surface = surf(gridX, gridY, gridP.');
    
    surface.EdgeColor = 'white';
    surface.EdgeAlpha = 0.25;
    
    % view(2)
    view([120 120 120])
    
    xlabel('x-dimension',     'Interpreter','latex')
    ylabel('y-dimension',     'Interpreter','latex')
    zlabel('pressure',     'Interpreter','latex')
    
    title(sprintf('t = %.2f s', t0), 'Interpreter','latex')
    t0 = t0 + dt;

    colormap viridis
    zlim([0 16])
    caxis([0 16])
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












