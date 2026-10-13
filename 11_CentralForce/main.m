clear; clc; close all;
setDefaults();

% Time
t0 = 0;
tf = 200;
tspan = [t0 tf];
speedupFactor = 50;

% System Params
m = 1; M = 1;
G = -1.5; k = 1;
r0 = [3; 0];
v0 = [0; 1];
z0 = [r0; v0];

fixedX0 = 2; %%% POINCARE SECTION

% Pack params struct
p.m = m; p.tspan = tspan; 
p.M = M; p.G = G; p.k = k;
p.z0 = z0;
p.fixedX0 = fixedX0;

% Root finding to find periodic solutions
guess_tf = 30;
[zstar, timePeriod] = doRootFinding(p, guess_tf);

zstar = [fixedX0; zstar];

if isempty(zstar)
    zstar = z0;
    % return;
end

z0 = zstar;
p.z0 = zstar;

% Find Jacobian

%%%%% Below is for finding jacobian of a fixed point
% h = 1e-3;
% jac = zeros(numel(zstar), numel(zstar));
% for i = 1:numel(zstar)
%     purturbation = zeros(numel(zstar), 1);
%     purturbation(i) = h;
%     zstar_plus = zstar + purturbation;
%     zstar_minus = zstar - purturbation;
% 
%     jac(:, i) = (centralForceODE(0, zstar_plus, p) - centralForceODE(0, zstar_minus, p)) / (2*h);
% end
% jac
% eigvals = eig(jac)

%%%%% Below is for finding jacobian of a period poincare map
T    = timePeriod;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
zAtT = @(z) deval(ode45(@(t,z) centralForceODE(t,z,p), [0 T], z, opts), T);

n   = numel(zstar);
jac = zeros(n);
h   = 1e-6;
for j = 1:n
    e = zeros(n,1);  e(j) = h;
    jac(:,j) = (zAtT(zstar + e) - zAtT(zstar - e)) / (2*h);
end
eigvals = eig(jac)
real(eigvals)


% Solve ODE
odeFunc = @(t, z) centralForceODE(t, z, p);
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-12);
sol = ode45(odeFunc, tspan, z0, odeOptions);

% Unpack sol
numPoints = 1e3 + 1;
tarray = linspace(t0, tf, numPoints);
zarray = deval(sol, tarray);
xarray = zarray(1, :);
yarray = zarray(2, :);
vxarray = zarray(3, :);
vyarray = zarray(4, :);

% Plot Solution
figure('Name', 'Trajectories')
tiledlayout(2, 4);

nexttile([1, 2]);
plot(xarray, yarray);
xlabel('x'); ylabel('y');
axis equal;

nexttile([1, 2]);
plot(vxarray, vyarray);
xlabel('v_x'); ylabel('v_y');
axis equal;

nexttile;
plot(tarray, xarray);
xlabel('t'); ylabel('x');

nexttile;
plot(tarray, yarray);
xlabel('t'); ylabel('y');

nexttile;
plot(tarray, vxarray);
xlabel('t'); ylabel('v_x');

nexttile;
plot(tarray, vyarray);
xlabel('t'); ylabel('v_y');

% Animation
figure('Name', 'Animation')
trajData = plot(xarray(1), yarray(1), 'r--'); hold on;
startPoint = plot(xarray(1), yarray(1), 'b.', 'MarkerSize', 20); hold on;
origin = plot(0, 0, 'k.', 'MarkerSize', 20);

axis equal;
axis([min(xarray) max(xarray) min(yarray) max(yarray)]);

tic;
t = 0;
while t <= tf 
    idx = find(tarray <= t, 1, "last");
    if ~isempty(idx)
        xData = xarray(1:idx);
        yData = yarray(1:idx);
        trajData.XData = xData;
        trajData.YData = yData;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end


function zdot = centralForceODE(t, z, p)
    m = p.m;
    M = p.M;
    G = p.G;
    k = p.k;

    x = z(1);
    y = z(2);
    vx = z(3);
    vy = z(4);
    i_hat = [1; 0];
    j_hat = [0; 1];

    r_vec = [x; y];
    r_mag = sqrt(x*x + y*y);
    r_hat = r_vec / r_mag;

    % F_r = G * m * M / r_mag^2; % F(r)
    % F_r = 1;
    F_r = G * m * M / r_mag^2 + k*r_mag; 
    centralForce = -F_r * r_hat;

    accel = centralForce ./ m;
    vx_dot = dot(accel, i_hat);
    vy_dot = dot(accel, j_hat);

    x_dot = vx;
    y_dot = vy;

    zdot = [x_dot; y_dot; vx_dot; vy_dot];
end

function [zstar, timePeriod] = doRootFinding(p, guess_tf)
    z0 = p.z0(2:end); tf = guess_tf;
    guess_vars = [z0; tf];
    rootFindFunc = @(guess_vars) findPeriodicSolutions(guess_vars, p);
    rootFindOptions = optimoptions('fsolve', 'Algorithm', ...
        'levenberg-marquardt', 'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 3000, 'FunctionTolerance', 1e-7);
    [z_vars,fval,exitflag,output] = fsolve(rootFindFunc, guess_vars, rootFindOptions)

    zstar = [];
    timePeriod = [];

    if exitflag > 0
        zstar = z_vars(1:end-1);
        timePeriod = z_vars(end);
    end
end

function error = findPeriodicSolutions(guess_vars, p)
    z0 = guess_vars(1:end-1);
    x0 = p.fixedX0;
    z0 = [x0; z0];
    tf = guess_vars(end);

    odeFunc = @(t, z) centralForceODE(t, z, p);
    sol = ode45(odeFunc, [p.tspan(1) tf], z0);
    error = z0 - deval(sol, tf);
end


function setDefaults()
    set(groot, 'defaultAxesXGrid', 'on', ...
               'defaultAxesYGrid', 'on', ...
               'defaultAxesZGrid', 'on')
end

