clear; close all; clc;

% System Params
p.m1 = 1;
p.m2 = 1;
p.k  = 100;
p.l0 = 1;
p.F0 = 1;

m1 = p.m1; m2 = p.m2; l0 = p.l0; F0 = p.F0;

% IC (both masses at rest, spring unstretched)
x1_0 = 0;  v1_0 = 0;
x2_0 = l0; v2_0 = 0;

z0 = [x1_0; v1_0; x2_0; v2_0];

% Time Params (Sim Params)
t0 = 0; tf = 10;
tspan = [t0 tf];
speedupFactor = 1;
p.tspan = tspan;

% Solve ODE
odefunc = @(t, z) twoMassesSpringRHS(t, z, p);
small = 1e-8;
odeOptions = odeset('RelTol', small, 'AbsTol', small);
sol = ode45(odefunc, tspan, z0, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray);

x1array = zarray(1, :);
v1array = zarray(2, :);
x2array = zarray(3, :);
v2array = zarray(4, :);
larray = x2array - x1array;


%%% Arrays for part B
% should your numerical solution always have the  
% property that F = m_tot a_G same as F0
xcomarray = (m1 * x1array + m2 * x2array) / (m1 + m2);
vcomarray = gradient(xcomarray, tarray);
acomarray = gradient(vcomarray, tarray);
F_G = (m1 + m2) * acomarray;

diffInF = F_G - F0;
diffInF = diffInF(3:end-2);
tarrayForDiffInF = tarray(3:end-2);

rmsError = sqrt(mean(diffInF .^ 2));
rmsErrorArray = rmsError * ones(size(tarrayForDiffInF));

%%% Arrays for part C
% For large enough K, the motion (displacment from initial position) is 
% close to CoM motion.
% Sweep k and record RMS deviation of each mass from CoM motion
kArray = logspace(0, 6, 20);
rmsDiffM1ComArray = zeros(size(kArray));
rmsDiffM2ComArray = zeros(size(kArray));

for i = 1:length(kArray)
    pk = p;
    pk.k = kArray(i);

    odefuncK = @(t, z) twoMassesSpringRHS(t, z, pk);
    solK = ode45(odefuncK, tspan, z0, odeOptions);
    zarrayK = deval(solK, tarray);

    x1arrayK = zarrayK(1, :);
    x2arrayK = zarrayK(3, :);
    xcomarrayK = (m1 * x1arrayK + m2 * x2arrayK) / (m1 + m2);

    diffM1ComK = (x1arrayK - x1_0) - (xcomarrayK - xcomarrayK(1));
    diffM2ComK = (x2arrayK - x2_0) - (xcomarrayK - xcomarrayK(1));
    rmsDiffM1ComArray(i) = sqrt(mean(diffM1ComK.^2));
    rmsDiffM2ComArray(i) = sqrt(mean(diffM2ComK.^2));
end

% Plot
figure('Name', 'Trajectory')
tiledlayout(2, 1);

nexttile;
plot(tarray, x1array, 'r-'); hold on;
plot(tarray, x2array, 'b-');
plot(tarray, xcomarray, 'k--');

legend('x_1', 'x_2', 'x_{COM}');
xlabel('t (s)');
ylabel('position (m)');
grid on;

nexttile;
plot(tarray, larray, 'k-'); hold on;
yline(l0, 'k--');

legend('l = x_2 - x_1', 'l_0');
xlabel('t (s)');
ylabel('spring length (m)');
grid on;

figure('Name', 'PartB')
plot(tarrayForDiffInF, diffInF, 'DisplayName', 'F_G - F_0'); hold on;
plot(tarrayForDiffInF, rmsErrorArray, '-.', 'DisplayName', "RMS Error = " + num2str(rmsError)); 
title('Comparison of $F_G = (m_1 + m_2) a_G$ and $F_0$', 'Interpreter', 'latex')
legend();

figure('Name', 'PartC')
loglog(kArray, rmsDiffM1ComArray, 'r.-', 'DisplayName', 'm_1'); hold on;
loglog(kArray, rmsDiffM2ComArray, 'b.-', 'DisplayName', 'm_2');
title('RMS deviation of mass displacement from CoM displacement')
xlabel('k (N/m)');
ylabel('RMS deviation (m)');
legend();
grid on;


% Animation
figure('Name', 'Animation');

springLine = plot([x1array(1) x2array(1)], [0 0], 'k-', 'LineWidth', 2); hold on;
mass1Data = plot(x1array(1), 0, 'r.', 'MarkerSize', 40 * sqrt(m1));
mass2Data = plot(x2array(1), 0, 'b.', 'MarkerSize', 40 * sqrt(m2));
comData = plot(xcomarray(1), 0, 'kx', 'MarkerSize', 12);

legend('Spring', 'm1', 'm2', 'COM');
xlabel('x (m)');
grid on;
axis([min([x1array x2array]) - l0, max([x1array x2array]) + l0, -1, 1]);
axis equal;

tic;
t = 0;
while t <= tf
    tnow_idx = find(tarray <= t, 1, "last");
    if ~isempty(tnow_idx)
        znow = zarray(:, tnow_idx);
        x1now = znow(1);
        x2now = znow(3);
        xcomnow = (m1 * x1now + m2 * x2now) / (m1 + m2);
        springLine.XData = [x1now x2now];
        mass1Data.XData = x1now;
        mass2Data.XData = x2now;
        comData.XData = xcomnow;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end
