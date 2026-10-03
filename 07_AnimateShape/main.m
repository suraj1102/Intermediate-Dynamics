clear; clc; close all;

% Coordinates of shape in object frame
shape.points = [
    [-1 ; 1], ... % top left
    [1 ; 1],  ... % top right
    [1 ;-1],  ... % bottom right
    [-1 ;-1], ... % bottom left
    ];

% Position in world coordinate frame
initialPosition = [0; 0];

% Transformation
T = [2; 2]; % translate by this much
theta = 20 * pi/180; % rotate by this much
R = [
        [cos(theta) -sin(theta)]; 
        [sin(theta) cos(theta)]
    ];

% Compute Transformation
transformedShape.points = [];
for i = 1:length(shape.points)
    transformedPoint = R * shape.points(:, i) + T;
    transformedShape.points = [transformedShape.points transformedPoint];
end


% Plot transformed shape
figure()
drawRecShape(shape); hold on;
drawRecShape(transformedShape, 'b', 'b');
axis equal;
grid on;
shg;


% Animation - linear interpolation between shapes
timeToAnimate = 4; % seconds to finish transformation.

deltaPoints = [];
for i = 1:length(shape.points)
    point = shape.points(:, i);
    transformed = transformedShape.points(:, i);
    delta = transformed - point;
    deltaPoints = [deltaPoints delta];
end


figure()
pointPlotDatas = [];
for i = 1:length(shape.points)
    pointPlotDatas = [pointPlotDatas plot(shape.points(1,:), shape.points(2,:), ...
        'r.', 'MarkerSize', 20);];
    hold on;
end

lineDatas = []';
for i = 1:length(shape.points)
    p1 = shape.points(:, i);
    p2 = shape.points(:, mod(i, length(shape.points)) + 1);
    lineDatas = [lineDatas line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'k')];
    hold on;
end

grid on;
axis equal;
axis([-5 5 -5 5]);

tic;
t_interpolate = 0;
while t_interpolate < timeToAnimate
    pointsAtT = shape.points + t_interpolate/timeToAnimate * deltaPoints;
    for i = 1:length(shape.points)
        pointPlotDatas(i).XData = pointsAtT(1, i);
        pointPlotDatas(i).YData = pointsAtT(2, i);
    end

    for i = 1:length(shape.points)
        lineDatas(i).XData = pointsAtT(:, i);
        lineDatas(i).YData = pointsAtT(:, mod(i, length(shape.points)) + 1);
        
    end


    t_interpolate = toc;
    drawnow limitrate
end




function drawRecShape(shape, pointColor, lineColor)
    %ASSUMES: shape.points arranged in circular manner.
    arguments
        shape
        pointColor = 'r'
        lineColor  = 'k'
    end
    
    n = size(shape.points, 2);
    
    % draw lines (including the closing edge back to the first point)
    for i = 1:n
        p1 = shape.points(:, i);
        p2 = shape.points(:, mod(i, n) + 1);
        line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', lineColor); hold on;
    end
    
    % draw points
    plot(shape.points(1,:), shape.points(2,:), '.', ...
         'Color', pointColor, 'MarkerSize', 20);
end