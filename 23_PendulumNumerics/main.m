clear; close all; clc;

% System Params
p.l = 1;
p.g = 10;
p.m = 1;

l = p.l; g = p.g; m = p.m;

% ICs
% theta measured from +x (horizontal), counter clockwise
% hanging straight down is theta = -pi/2
%%% LMB States
theta0 = deg2rad(30);
thetadot0 = 0;

z0LMB = [theta0; thetadot0];

%%% DAE States
% r = l [cos(theta); sin(theta)]
% v = rdot = l thetadot [-sin(theta); cos(theta)]
x0 = l * cos(theta0);
y0 = l * sin(theta0);
vx0 = -l * thetadot0 * sin(theta0);
vy0 = l * thetadot0 * cos(theta0);

z0DAE = [x0; y0; vx0; vy0];

% Time Params
t0 = 0;
tf = 100;
tspan = [t0 tf]; p.tspan = tspan;
speedupFactor = 1;

% Solve ODE
small = 1e-5;
odeOptions = odeset('RelTol', small, 'AbsTol', small);

odefuncLMB = @(t, z) pendulumLMBrhs(t, z, p);
solLMB = ode45(odefuncLMB, tspan, z0LMB, odeOptions);

odefuncDAE = @(t, z) pendulumDAErhs(t, z, p);
solDAE = ode45(odefuncDAE, tspan, z0DAE, odeOptions);

% Unpack ODE sol
tarray = linspace(t0, tf, 1e4+1);
zarrayLMB = deval(solLMB, tarray);
zarrayDAE = deval(solDAE, tarray);

thetaarray = zarrayLMB(1, :);
xarrayLMB = l * cos(thetaarray);
yarrayLMB = l * sin(thetaarray);

xarrayDAE = zarrayDAE(1, :);
yarrayDAE = zarrayDAE(2, :);
vxarrayDAE = zarrayDAE(3, :);
vyarrayDAE = zarrayDAE(4, :);

%%% Arrays for LMB vs DAE comparison
absErrorX = abs(xarrayDAE - xarrayLMB);
absErrorY = abs(yarrayDAE - yarrayLMB);
rmsErrorX = sqrt(mean(absErrorX.^2));
rmsErrorY = sqrt(mean(absErrorY.^2));

%%% Arrays for drift in constraint
% length constraint: sqrt(x^2 + y^2) = l
errorInLengthArray = sqrt(xarrayDAE.^2 + yarrayDAE.^2) - l;
% velocity constraint: r . v = 0
errorInRdotVArray = xarrayDAE .* vxarrayDAE + yarrayDAE .* vyarrayDAE;

% Plot
figure('Name', 'Trajectory')
tiledlayout(2, 1);

nexttile;
plot(tarray, xarrayLMB, 'r-'); hold on;
plot(tarray, xarrayDAE, 'b--');
legend('x LMB', 'x DAE');
xlabel('t (s)');
ylabel('x (m)');
grid on;

nexttile;
plot(tarray, yarrayLMB, 'r-'); hold on;
plot(tarray, yarrayDAE, 'b--');
legend('y LMB', 'y DAE');
xlabel('t (s)');
ylabel('y (m)');
grid on;

figure('Name', 'LMBvsDAE')
semilogy(tarray, absErrorX, 'r-', 'DisplayName', "|x_{DAE} - x_{LMB}|, RMS = " + num2str(rmsErrorX)); hold on;
semilogy(tarray, absErrorY, 'b-', 'DisplayName', "|y_{DAE} - y_{LMB}|, RMS = " + num2str(rmsErrorY));
title('Difference between DAE and LMB solutions')
xlabel('t (s)');
grid on;
legend();

figure('Name', 'ConstraintDrift')
tiledlayout(2, 1);

nexttile;
plot(tarray, errorInLengthArray, 'k-');
title('Error in length $= \sqrt{x^2 + y^2} - l$', 'Interpreter', 'latex')
xlabel('t (s)');
grid on;

nexttile;
plot(tarray, errorInRdotVArray, 'k-');
title('Error in velocity constraint $= \vec{r} \cdot \vec{v}$', 'Interpreter', 'latex')
xlabel('t (s)');
grid on;

% Save Figures
exportgraphics(findobj('Name', 'Trajectory'), 'imgs/trajectory.png', 'Resolution', 300);
exportgraphics(findobj('Name', 'LMBvsDAE'), 'imgs/LMBvsDAE.png', 'Resolution', 300);
exportgraphics(findobj('Name', 'ConstraintDrift'), 'imgs/constraintDrift.png', 'Resolution', 300);

% Animation
figure('Name', 'Animation');

rodLMB = plot([0 xarrayLMB(1)], [0 yarrayLMB(1)], 'r-', 'LineWidth', 2); hold on;
rodDAE = plot([0 xarrayDAE(1)], [0 yarrayDAE(1)], 'b--', 'LineWidth', 2);
bobLMB = plot(xarrayLMB(1), yarrayLMB(1), 'r.', 'MarkerSize', 40);
bobDAE = plot(xarrayDAE(1), yarrayDAE(1), 'bo', 'MarkerSize', 12);
plot(0, 0, 'k.', 'MarkerSize', 20);

legend('Rod LMB', 'Rod DAE', 'LMB', 'DAE', 'Pivot');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;
axis([-1.2 * l, 1.2 * l, -1.2 * l, 1.2 * l]);

tic;
t = 0;
while t <= tf
    tnow_idx = find(tarray <= t, 1, "last");
    if ~isempty(tnow_idx)
        xLMBnow = xarrayLMB(tnow_idx);
        yLMBnow = yarrayLMB(tnow_idx);
        xDAEnow = xarrayDAE(tnow_idx);
        yDAEnow = yarrayDAE(tnow_idx);
        rodLMB.XData = [0 xLMBnow];
        rodLMB.YData = [0 yLMBnow];
        rodDAE.XData = [0 xDAEnow];
        rodDAE.YData = [0 yDAEnow];
        bobLMB.XData = xLMBnow;
        bobLMB.YData = yLMBnow;
        bobDAE.XData = xDAEnow;
        bobDAE.YData = yDAEnow;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end

% Save final frame of animation (t = tf)
rodLMB.XData = [0 xarrayLMB(end)];
rodLMB.YData = [0 yarrayLMB(end)];
rodDAE.XData = [0 xarrayDAE(end)];
rodDAE.YData = [0 yarrayDAE(end)];
bobLMB.XData = xarrayLMB(end);
bobLMB.YData = yarrayLMB(end);
bobDAE.XData = xarrayDAE(end);
bobDAE.YData = yarrayDAE(end);
title("t = " + num2str(tarray(end)) + " s, length error = " + num2str(errorInLengthArray(end)));
drawnow;
exportgraphics(findobj('Name', 'Animation'), 'imgs/animationFinalFrame.png', 'Resolution', 300);
