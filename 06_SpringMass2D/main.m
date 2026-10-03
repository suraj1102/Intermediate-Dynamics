close all; clear; clc;

% System Params
p.m = 1;
p.g = 10;
p.k = 10;
p.l0 = 0.0;

[m, g, k, l0] = unpack_p(p);

% IC
r0 = [1; 0];
v0 = [0; 0];
z0 = [r0; v0];

% Time Params (Sim Params)
t0 = 0; tf = 3;
tspan = [t0 tf];
speedupFactor = 1;
p.tspan = tspan;

% Solve ODE
odefunc = @(t, z) springMassRHS(t, z, p);
odeOptions = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
sol = ode45(odefunc, tspan, z0, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray); % 4 x length(tarray)

xarray = zarray(1, :);
yarray = zarray(2, :);
vxarray = zarray(3, :);
vyarray = zarray(4, :);
origin = [0; 0];

% Plot
figure('Name', 'Trajectory')

plot(xarray, yarray, 'r-'); hold on;
plot(xarray(1), yarray(1), 'b.', 'MarkerSize', 20);
plot(xarray(end), yarray(end), 'g.', 'MarkerSize', 20);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);
drawCircle([0;0], l0);
plot([origin(1) xarray(end)], [origin(2) yarray(end)], 'k-', 'LineWidth', 1.5);

legend('Trajectory', 'Start', 'End', 'Origin', 'l0', 'Spring');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

% Animation
figure('Name', 'Animation');

trajData = plot(xarray(1), yarray(1), 'r-'); hold on;
massData = plot(xarray(1), yarray(1), 'b.', 'MarkerSize', 20);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);
springLine = plot([origin(1) xarray(1)], [origin(2) yarray(1)], 'k-', 'LineWidth', 1.5);
drawCircle([0;0], l0);

axis equal;
grid on;
axis([min([xarray 0]) max([xarray 0]) min([yarray 0]) max([yarray 0])]);

tic;
t = 0;
while t <= tf
    tnow_idx = find(tarray <= t, 1, "last");
    if ~isempty(tnow_idx)
        znow = zarray(:, tnow_idx);
        xnow = znow(1);
        ynow = znow(2);
        trajData.XData = [trajData.XData xnow];
        trajData.YData = [trajData.YData ynow];
        massData.XData = xnow;
        massData.YData = ynow;
        springLine.XData = [origin(1) xnow];
        springLine.YData = [origin(2) ynow];
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end

