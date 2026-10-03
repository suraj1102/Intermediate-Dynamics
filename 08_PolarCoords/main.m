clear; clc; close all;

% Time Params
t0 = 0;
tf = 10;
tspan = [t0 tf];
numPoints = 1e3+1;
speedupFactor = 5;

% System Params
p = struct();

% System IC

z0 = [r0; theta0; rdot0; thetadot0];

% Solve ODE
odeFunc = @(t, z) myRHS(t, z, p);
sol = ode45(odeFunc, tspan, z0);

tarray = linspace(t0, tf, numPoints);
zarray = deval(sol, tarray);            % 4 x numPoints
r   = zarray(1,:);  th  = zarray(2,:);
rd  = zarray(3,:);  thd = zarray(4,:);
cart = polarToCartesian(zarray);        % 4 x numPoints
x  = cart(1,:);  y  = cart(2,:);
vx = cart(3,:);  vy = cart(4,:);

% % % Plot Trajectory
figure('Position', [100 100 1400 800]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% x vs y, vx vs vy
axXY = nexttile([1 2]);
plot(x, y, 'b'); hold on; grid on; axis equal;
xlabel('x'); ylabel('y'); title('Position');
posMarker = plot(x(1), y(1), 'r.', 'MarkerSize', 25);

axV = nexttile([1 2]);
plot(vx, vy, 'b'); hold on; grid on; axis equal;
xlabel('v_x'); ylabel('v_y'); title('Velocity');
velMarker = plot(vx(1), vy(1), 'r.', 'MarkerSize', 25);

% states vs time
data   = {x, y, vx, vy, r, th, rd, thd};
labels = {'x', 'y', 'v_x', 'v_y', 'r', '\theta', 'dr/dt', 'd\theta/dt'};
cursors = gobjects(1, numel(data));
for k = 1:numel(data)
    nexttile;
    plot(tarray, data{k}, 'b'); hold on; grid on;
    xlabel('t'); ylabel(labels{k});
    cursors(k) = xline(t0, 'r');        % moving time cursor
end

%% Animate
tic;
t = 0;
while t < tf
    [~, idx] = min(abs(tarray - t));    % closest sample to current time
    posMarker.XData = x(idx);   posMarker.YData = y(idx);
    velMarker.XData = vx(idx);  velMarker.YData = vy(idx);
    set(cursors, 'Value', tarray(idx));
    title(axXY, sprintf('Position   (t = %.2f s)', tarray(idx)));
    drawnow limitrate
    t = toc * speedupFactor;
end

%% Straightness study
%% Error vs analytical solution
% --- ode45: sweep tolerances
tols = 10.^(-3:-1:-12);
err45 = zeros(size(tols));
for k = 1:numel(tols)
    opts = odeset('RelTol', tols(k), 'AbsTol', tols(k) * 1e-2);
    [t, z] = ode45(odeFunc, tspan, z0, opts);
    err45(k) = trueError(t, z, z0);
end
% --- Forward Euler: sweep step size
hs = 10.^(-1:-0.5:-5);
errEuler = zeros(size(hs));
for k = 1:numel(hs)
    [t, z] = eulerSolver(odeFunc, tspan, z0, odeset('MaxStep', hs(k)));
    errEuler(k) = trueError(t, z, z0);
end
% --- Plots
figure('Position', [100 100 1400 450]);
tiledlayout(1, 3, 'TileSpacing', 'compact');
nexttile;
loglog(hs, errEuler, 'o-'); hold on; grid on;
loglog(hs, errEuler(end) * hs / hs(end), 'k--');   % slope-1 reference
xlabel('step size h'); ylabel('max distance from exact solution');
title('Forward Euler'); legend('error', 'O(h) reference', 'Location', 'northwest');
nexttile;
loglog(tols, err45, 'o-'); hold on; grid on;
loglog(tols, tols, 'k--');                          % error = tol reference
xlabel('RelTol'); ylabel('max distance from exact solution');
title('ode45'); legend('error', 'error = tol', 'Location', 'northwest');
nexttile; hold on; grid on; axis equal;
for h = [1e-1 1e-2 1e-3]
    [t, z] = eulerSolver(odeFunc, tspan, z0, odeset('MaxStep', h));
    c = polarToCartesian(z.');
    plot(c(1,:), c(2,:), 'DisplayName', sprintf('Euler h = %g', h));
end
[t, z] = ode45(odeFunc, tspan, z0, odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
c = polarToCartesian(z.');
plot(c(1,:), c(2,:), 'DisplayName', 'ode45 tol 1e-10');
c0 = polarToCartesian(z0);
tEx = [t0 tf];
plot(c0(1) + c0(3) * tEx, c0(2) + c0(4) * tEx, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'exact');
xlabel('x'); ylabel('y'); title('Paths'); legend('Location', 'best');

%% Function Defs
function zdot = myRHS(t, z, p)
    r = z(1);
    theta = z(2);
    rdot = z(3);
    thetadot = z(4);

    rddot = r * thetadot^2;
    thetaddot = -2*rdot*thetadot/r;

    zdot = [rdot; thetadot; rddot; thetaddot];
end

function cart = polarToCartesian(z)
    % z: 4 x N polar state [r; theta; rdot; thetadot] -> 4 x N [x; y; vx; vy]
    r = z(1,:); th = z(2,:); rd = z(3,:); thd = z(4,:);
    x  = r .* cos(th);
    y  = r .* sin(th);
    vx = rd .* cos(th) - r .* thd .* sin(th);
    vy = rd .* sin(th) + r .* thd .* cos(th);
    cart = [x; y; vx; vy];
end


function err = trueError(t, z, z0)
% Max distance between the numerical path and the exact free-particle
% solution x(t) = x0 + v0*(t - t0), evaluated at the solver's own times.
% t: N x 1, z: N x 4 polar states (same layout as ode45 output).
    c0 = polarToCartesian(z0);                       % [x0; y0; vx0; vy0]
    cart = polarToCartesian(z.');                    % 4 x N
    exact = c0(1:2) + c0(3:4) .* (t(:).' - t(1));    % 2 x N
    err = max(vecnorm(cart(1:2, :) - exact));
end


function [tout, yout] = eulerSolver(odefun, tspan, y0, options)
% Fixed-step forward Euler with the same call signature as ode45.
% Step size comes from odeset('MaxStep', h); default is (tf - t0)/1000.
    t0 = tspan(1);  tf = tspan(end);
    h = (tf - t0) / 1000;
    if nargin >= 4 && ~isempty(odeget(options, 'MaxStep'))
        h = odeget(options, 'MaxStep');
    end
    N = ceil((tf - t0) / h * (1 - 1e-10));   % relative tolerance
    tout = zeros(N + 1, 1);
    yout = zeros(N + 1, numel(y0));
    tout(1) = t0;  yout(1, :) = y0(:).';
    y = y0(:);
    for n = 1:N
        tNext = min(t0 + n * h, tf);
        if tNext <= tout(n)                  % already at tf: trim and stop
            tout = tout(1:n);  yout = yout(1:n, :);
            break
        end
        y = y + (tNext - tout(n)) * odefun(tout(n), y);
        tout(n + 1) = tNext;
        yout(n + 1, :) = y.';
    end
    tout(end) = tf;                          % land exactly on tf
end