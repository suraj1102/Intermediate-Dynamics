clear; clc; close all;

%% Parameters

% Time params
t0 = 0;
tf = 10;
tspan = [t0 tf];
h = 1e-3;

numPoints = (tf - t0) / h + 1;
speedupFactor = 0.5;

% System params
m = 10; % kg
k = 100; % N/m
l0 = 2; % m
b = 10; % Ns/m

x0 = 1; % m
vx0 = 1; % m/s
z0 = [x0 vx0]; % Initial Conditions

% Params struct
p.m = m;
p.k = k;
p.l0 = l0;
p.b = b;
p.tspan = tspan;

%% ODE Solving
odefunc = @(t, z) myRHS(t, z, p);
[tarray, zarray] = eulerSolver(odefunc, tspan, z0, h);

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
plot(l0, linspace(-1, 1, 1000), 'b--');

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

%% Functions

function [t, z] = eulerSolver(odefunc, tspan, z0, h)
    t = tspan(1) : h : tspan(2);
    t = t'; % column vector
    z = zeros(size(t, 1), length(z0)); % N, state_dim

    z(1, :) = z0;

    for i = 1 : size(z, 1) - 1
        t_now = t(i);
        z_now = z(i, :);
        z_next = z_now + h * odefunc(t_now, z_now);

        z(i+1, :) = z_next;
    end
end



function zdot = myRHS(t, z, p)
    l0 = p.l0;
    b = p.b;
    m = p.m;
    k = p.k;

    x = z(1);
    vx = z(2);

    xdot = vx;
    vxdot = -k/m * (x - l0) - b/m * vx;

    zdot = [xdot, vxdot];
end

