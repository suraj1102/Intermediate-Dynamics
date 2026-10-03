clear; clc; close all;

%% Parameters

% Time params
t0 = 0;
tf = 10;
tspan = [t0 tf];
h = 1e-3;

numPoints = (tf - t0) / h + 1;
speedupFactor = 1;

% System params
m = 10; % kg
k = 100; % N/m
l0 = 2; % m
b = 10; % Ns/m

x0 = 1; % m
vx0 = 1; % m/s
z0 = [x0 vx0]; % Initial Conditions

% Params struct
p.m = 10;
p.tspan = tspan;

%% ODE Solving
tarray = linspace(t0, tf, numPoints)';
zarray = zeros(length(tarray), length(z0));
zarray(1, :) = z0;

for i = 1 : numPoints - 1
    t_now = tarray(i);
    x_now = zarray(i, 1);
    vx_now = zarray(i, 2);
    
    xdot_now = vx_now;
    x_next = x_now + h * xdot_now;

    vxdot_now = -k/m * (x_now - l0) -b/m * vx_now;
    vx_next = vx_now + h * vxdot_now;

    z_next = [x_next vx_next];
    zarray(i+1, :) = z_next;
end

%% Plotting
xarray = zarray(:, 1);
vxarray = zarray(:, 2);

figure()
tiledlayout(2, 1);

nexttile;
plot(tarray, xarray, 'DisplayName', 'x(t)')
grid on;
xlabel("t")
ylabel("x")

nexttile;
plot(tarray, vxarray, 'DisplayName', 'vx(t)')
grid on;
xlabel("t")
ylabel("vx")

%% Animation

figure('Name', 'Animaton')
trajData = plot(xarray(1), 0, 'r.', 'MarkerSize', 20);
hold on;
plot(l0 * ones(1, 1000), linspace(-1, 1, 1000), 'b--');

grid on;
axis equal;
axis([min(xarray), max(xarray), -0.3, 0.3])


tic;
t = 0;
while t < tf

    idx = find(tarray <= t, 1, 'last');
    
    if ~isempty(idx)
        trajData.XData = xarray(idx);
    end

    drawnow limitrate;

    t = toc * speedupFactor;
end