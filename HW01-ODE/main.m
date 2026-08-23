clear; close all; clc;

p = struct();

t0 = 0;
tf = 15;
tspan = [t0 tf];
h = 0.01;

z0 = [0, 0, 0.5, -0.5];


odefunc = @(t, z) myRHS(t, z, p);
[t_all, z_all] = eulerSolver(odefunc, tspan, z0, h);


x_all  = z_all(:, 1);
y_all  = z_all(:, 2);
vx_all = z_all(:, 3);
vy_all = z_all(:, 4);


figure('Name', "States")
tiledlayout(2, 2);


nexttile;
plot(t_all, x_all);
title("x(t)")
grid on;

nexttile;
plot(t_all, y_all);
title("y(t)")
grid on;

nexttile;
plot(t_all, vx_all);
title("vx(t)")
grid on;

nexttile;
plot(t_all, vy_all);
title("vy(t)")
grid on;

figure('Name', "Phase Plot")
plot(x_all(1), y_all(1), 'DisplayName', 'Start', 'Marker', '.', ...
    'MarkerSize', 20); hold on;
plot(x_all(end), y_all(end), 'DisplayName', 'End', 'Marker', '.', ...
    'MarkerSize', 20); hold on;
plot(x_all, y_all, 'DisplayName', 'Trajectory')
legend()
axis equal
grid on;
shg


figure("Name", "Animation")


trajStatic = plot(x_all, y_all, 'r--', 'HandleVisibility', 'off');
hold on;
trajData = plot(x_all(1), y_all(1), 'b-', 'HandleVisibility','off');
hold on;
currentPoint = plot(x_all(1), y_all(1), 'g.', 'MarkerSize', 20, 'HandleVisibility','off');


xlabel("x");
ylabel("y");

axis equal;
axis([min(x_all) max(x_all) min(y_all) max(y_all)]);
grid on;
legend off;



tic;
tSim = 0;
speedupFactor = 2;

while tSim <= tf
    % Find idx of tSim in tarray (t)
    idx = find(t_all <= tSim, 1, 'last');

    if ~isempty(idx)
        trajData.XData = x_all(1:idx);
        trajData.YData = y_all(1:idx);

        currentPoint.XData = x_all(idx);
        currentPoint.YData = y_all(idx);
    end

    drawnow limitrate

    tSim = toc * speedupFactor;
end





