clear; close all; clc;

% System Params
p.m1 = 1;
p.m2 = 1;
p.l1 = 1;
p.l2 = 1;
p.g = 10;

m1 = p.m1; m2 = p.m2; l1 = p.l1; l2 = p.l2; g = p.g;

% ICs
theta1 = pi/2;
theta2 = 0;
theta1dot = 0;
theta2dot = 1;

rE0 = l1 * [cos(theta1); sin(theta1)];
vE0 = l1 * theta1dot * [-sin(theta1); cos(theta1)];
rP0 = rE0 + l2 * [cos(theta2); sin(theta2)];
vP0 = vE0 + l2 * theta2dot * [-sin(theta2); cos(theta2)];

z0 = [rE0; vE0; rP0; vP0];

% Time Params
t0 = 0;
tf = 10;
tspan = [t0 tf]; p.tspan = tspan;
speedupFactor = 1;

% Solve ODE
odefunc = @(t, z) doublePendDAErhs(t, z, p);
small = 1e-5;
odeOptions = odeset('RelTol', small, 'AbsTol', small);
sol = ode45(odefunc, tspan, z0, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e4+1);
zarray = deval(sol, tarray);

xEarray = zarray(1, :);
yEarray = zarray(2, :);
xParray = zarray(5, :);
yParray = zarray(6, :);
origin = [0; 0];

% Plot
figure('Name', 'Trajectory')
plot(xParray, yParray, 'b-'); hold on;
plot(xEarray, yEarray, 'r-');
plot(xParray(1), yParray(1), 'b.', 'MarkerSize', 20);
plot(xEarray(1), yEarray(1), 'r.', 'MarkerSize', 20);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);

legend('P (forearm end)', 'E (elbow)', 'P Start', 'E Start', 'O');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

% Animation
figure('Name', 'Animation');

trajPData = plot(xParray(1), yParray(1), 'b-'); hold on;
armLine = plot([origin(1) xEarray(1) xParray(1)], [origin(2) yEarray(1) yParray(1)], 'k-', 'LineWidth', 2);
massEData = plot(xEarray(1), yEarray(1), 'r.', 'MarkerSize', 30);
massPData = plot(xParray(1), yParray(1), 'b.', 'MarkerSize', 30);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);

legend('P Trajectory', 'Arms', 'E', 'P', 'O');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;
axis([-1.1, 1.1, -1.1, 1.1] * (l1 + l2));

tic;
t = 0;
while t <= tf
    tnow_idx = find(tarray <= t, 1, "last");
    if ~isempty(tnow_idx)
        xEnow = xEarray(tnow_idx);
        yEnow = yEarray(tnow_idx);
        xPnow = xParray(tnow_idx);
        yPnow = yParray(tnow_idx);
        trajPData.XData = xParray(1:tnow_idx);
        trajPData.YData = yParray(1:tnow_idx);
        armLine.XData = [origin(1) xEnow xPnow];
        armLine.YData = [origin(2) yEnow yPnow];
        massEData.XData = xEnow;
        massEData.YData = yEnow;
        massPData.XData = xPnow;
        massPData.YData = yPnow;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end
