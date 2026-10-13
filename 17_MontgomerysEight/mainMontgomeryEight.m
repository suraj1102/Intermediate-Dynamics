close all; clear; clc;

% System Params
p.G = 1;
p.m1 = 1;
p.m2 = 1;
p.m3 = 1;

% IC

% GIVEN IC
perturb = 0.0;
r1_0 = [-0.97000436; 0.24308753];
r2_0 = -r1_0;
r3_0 = [0; 0];
v3_0 = [0.93240737; 0.86473146] + [perturb; 0];
v1_0 = -v3_0 / 2;
v2_0 = -v3_0 / 2;

% OVERWRITE IC
% r1_0 = [1; 0];
% r2_0 = [0; 1];
% r3_0 = [0; 0];
% v1_0 = [0; 1];
% v2_0 = [-1; 0];
% v3_0 = [0; 0];

z0 = [r1_0; r2_0; r3_0; v1_0; v2_0; v3_0];

% Time Params (Sim Params)
t0 = 0; tf = 10;
tspan = [t0 tf];
speedupFactor = 1;
p.tspan = tspan;

% Solve ODE
odefunc = @(t, z) myRHS(t, z, p);
odeOptions = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
sol = ode45(odefunc, tspan, z0, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray);

x1array = zarray(1, :);
y1array = zarray(2, :);
x2array = zarray(3, :);
y2array = zarray(4, :);
x3array = zarray(5, :);
y3array = zarray(6, :);
origin = [0; 0];

% Plot
figure('Name', 'Trajectory')
tiledlayout(1, 2);

nexttile;
plot(x1array, y1array, 'r-'); hold on;
plot(x2array, y2array, 'b-');
plot(x3array, y3array, 'g-');
plot(x1array(1), y1array(1), 'r.', 'MarkerSize', 20);
plot(x2array(1), y2array(1), 'b.', 'MarkerSize', 20);
plot(x3array(1), y3array(1), 'g.', 'MarkerSize', 20);

legend('m1', 'm2', 'm3', 'm1 Start', 'm2 Start', 'm3 Start');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

nexttile;
plot(tarray, x1array, 'r-'); hold on;
plot(tarray, x2array, 'b-');
plot(tarray, x3array, 'g-');

legend('x_1', 'x_2', 'x_3');
xlabel('t (s)');
ylabel('x (m)');
grid on;

% Animation
figure('Name', 'Animation');

traj1Data = plot(x1array(1), y1array(1), 'r-'); hold on;
traj2Data = plot(x2array(1), y2array(1), 'b-');
traj3Data = plot(x3array(1), y3array(1), 'g-');
mass1Data = plot(x1array(1), y1array(1), 'r.', 'MarkerSize', 20);
mass2Data = plot(x2array(1), y2array(1), 'b.', 'MarkerSize', 20);
mass3Data = plot(x3array(1), y3array(1), 'g.', 'MarkerSize', 20);
plot(origin(1), origin(2), 'k.', 'MarkerSize', 10);

axis equal;
grid on;
axis([min([x1array x2array x3array]) max([x1array x2array x3array]) min([y1array y2array y3array]) max([y1array y2array y3array])]);

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
        x3now = znow(5);
        y3now = znow(6);
        traj1Data.XData = [traj1Data.XData x1now];
        traj1Data.YData = [traj1Data.YData y1now];
        traj2Data.XData = [traj2Data.XData x2now];
        traj2Data.YData = [traj2Data.YData y2now];
        traj3Data.XData = [traj3Data.XData x3now];
        traj3Data.YData = [traj3Data.YData y3now];
        mass1Data.XData = x1now;
        mass1Data.YData = y1now;
        mass2Data.XData = x2now;
        mass2Data.YData = y2now;
        mass3Data.XData = x3now;
        mass3Data.YData = y3now;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end
