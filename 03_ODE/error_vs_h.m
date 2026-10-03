clear; clc; close all;

%% Parameters

% Time params
t0 = 0;
tf = 10;
tspan = [t0 tf];
h = 10.^(-8:1:1);

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

z_final = zeros(length(h), length(z0)); 

t_all = cell(length(h), 1);
z_all = cell(length(h), 1);

for i = 1:length(h)
    [tarray, zarray] = eulerSolver(odefunc, tspan, z0, h(i));
    
    t_all{i} = tarray;
    z_all{i} = zarray;
    
    z_final(i, :) = zarray(end, :);
end

%% Plotting Error vs. h
% Treat the finest grid (h = 10^-5, which is index 1) as the "exact" solution
x_exact = z_final(1, 1); 

% Calculate the absolute error of final position for all other step sizes
errors = abs(z_final(:, 1) - x_exact);

% We drop the first point for plotting because log10(0) is undefined
h_plot = h(2:end);
err_plot = errors(2:end);

figure('Name', 'Error Analysis');
loglog(h_plot, err_plot, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Euler Error');
hold on;

% Add a reference line for O(h) convergence (Slope = 1)
% We anchor it to the error at h = 0.1 for visual alignment
anchor_idx = find(h_plot == 0.1, 1);
if isempty(anchor_idx), anchor_idx = length(h_plot); end
ref_line = err_plot(anchor_idx) * (h_plot / h_plot(anchor_idx));

loglog(h_plot, ref_line, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h) Slope = 1');

grid on;
xlabel('Step Size (h)');
ylabel('Absolute Error at t = t_f');
title('Global Error of Euler Method');
legend('Location', 'northwest');



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

