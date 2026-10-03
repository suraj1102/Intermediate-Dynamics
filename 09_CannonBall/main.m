close all; clear; clc;

% System Params
p.m = 1;
p.g = 10;
p.c = 10; % Drag

[m, g, c] = unpack_p(p);

% IC
theta0 = deg2rad(45);
v0_mag = 10;
r0 = [0; 0];
v0 = [v0_mag * cos(theta0); v0_mag * sin(theta0)];
z0 = [r0; v0];

% Time Params (Sim Params)
t0 = 0; tf = 1000;
tspan = [t0 tf];
speedupFactor = 1;
p.tspan = tspan;

% Solve ODE
odefunc = @(t, z) cannonBallRHS(t, z, p);
eventfunc = @(t, z) landEvent(t, z);

odeOptions = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Events', eventfunc);
sol = ode45(odefunc, tspan, z0, odeOptions);

tf = sol.x(end);
tpsan = [t0 tf];
p.tspan = tspan;


% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray); % 4 x length(tarray)

xarray = zarray(1, :);
yarray = zarray(2, :);
vxarray = zarray(3, :);
vyarray = zarray(4, :);

% Plot
figure('Name', 'Trajectory')

plot(xarray, yarray, 'r-'); hold on;

grid on;
legend('Trajectory');

axis equal;

% Animation
figure('Name', 'Animation');

trajData = plot(xarray(1), yarray(1), 'r-'); hold on;
massData = plot(xarray(1), yarray(1), 'b.', 'MarkerSize', 20);

grid on;
axis equal;
axis([min(xarray) max(xarray) min(yarray) max(yarray)]);

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
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end

function [value,isterminal,direction] = landEvent(t,z)
    value = z(2); % value we want to be zero, y position here
    isterminal = 1; % yes, terminate ode integration
    direction = -1; % when y coord is decreasing and reaches zero, trigger event.
end

