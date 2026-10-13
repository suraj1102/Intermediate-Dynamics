close all; clear; clc;

% System Params
p.G = 100;
p.m1 = 10000;
p.m2 = 10;

G = p.G; m1 = p.m1; m2 = p.m2;

% IC
r1_0 = [0; 0];
r2_0 = [1; 0];


d = norm(r2_0 - r1_0);
omega = sqrt(G * (m1 + m2) / d^3);
T = 2*pi / omega;

r_com = (m1 * r1_0 + m2 * r2_0) / (m1 + m2);
v1_0 = [0; omega * (r1_0(1) - r_com(1))];
v2_0 = [0; omega * (r2_0(1) - r_com(1))];

z0 = [r1_0; r2_0; v1_0; v2_0];

% Time Params (Sim Params)
t0 = 0; tf = 3 * T;
tspan = [t0 tf];
speedupFactor = 1;
p.tspan = tspan;

% Solve ODE
odefunc = @(t, z) inverseGravityRHS(t, z, p);
odeOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
sol = ode45(odefunc, tspan, z0, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray);

x1array = zarray(1, :);
y1array = zarray(2, :);
x2array = zarray(3, :);
y2array = zarray(4, :);
darray = sqrt((x2array - x1array).^2 + (y2array - y1array).^2);
origin = [0; 0];

% Plot
figure('Name', 'Trajectory')
tiledlayout(1, 2);

nexttile;
plot(x1array, y1array, 'r-'); hold on;
plot(x2array, y2array, 'b-');
plot(x1array(1), y1array(1), 'r.', 'MarkerSize', 20);
plot(x2array(1), y2array(1), 'b.', 'MarkerSize', 20);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);

legend('m1', 'm2', 'm1 Start', 'm2 Start', 'COM');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

nexttile;
plot(tarray, x1array, 'r-'); hold on;
plot(tarray, x2array, 'b-');
plot(tarray, darray, 'k-');

legend('x_1', 'x_2', 'd');
xlabel('t (s)');
ylabel('position (m)');
grid on;

% Animation
figure('Name', 'Animation');

traj1Data = plot(x1array(1), y1array(1), 'r-'); hold on;
traj2Data = plot(x2array(1), y2array(1), 'b-');
mass1Data = plot(x1array(1), y1array(1), 'r.', 'MarkerSize', 20 * sqrt(m1));
mass2Data = plot(x2array(1), y2array(1), 'b.', 'MarkerSize', 20 * sqrt(m2));
plot(origin(1), origin(2), 'k.', 'MarkerSize', 20);
gravityLine = plot([x1array(1) x2array(1)], [y1array(1) y2array(1)], 'k--');

axis equal;
grid on;
axis([min([x1array x2array]) max([x1array x2array]) min([y1array y2array]) max([y1array y2array])]);

tic;
t = 0;
while t <= tf
    tnow_idx = find(tarray <= t, 1, "last");
    if ~isempty(tnow_idx)
        znow = zarray(:, tnow_idx);
        x1now = znow(1);
        y1now = znow(2);
        x2now = znow(3);
        y2now = znow(4);
        traj1Data.XData = [traj1Data.XData x1now];
        traj1Data.YData = [traj1Data.YData y1now];
        traj2Data.XData = [traj2Data.XData x2now];
        traj2Data.YData = [traj2Data.YData y2now];
        mass1Data.XData = x1now;
        mass1Data.YData = y1now;
        mass2Data.XData = x2now;
        mass2Data.YData = y2now;
        gravityLine.XData = [x1now x2now];
        gravityLine.YData = [y1now y2now];
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end
