clear; close all; clc;

% System Params
p.m1 = 1;
p.m2 = 1;
p.k  = 1e5;
p.l0 = 1;
p.d = p.l0;
p.F0 = 1;
% p.F0_func = @(t) constantForceFunction(t, p);
p.F0_func = @(t) sinusoidalForceFunction(t, p);

m1 = p.m1; m2 = p.m2; l0 = p.l0; F0 = p.F0; d = p.d;

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
odefunc = @(t, z) twoMassesConstraintDAE(t, z, p);
small = 1e-8;
odeOptions = odeset('RelTol', small, 'AbsTol', small);
ode45(odefunc, tspan, z0, odeOptions); % warm-up so timing excludes first-call overhead
tic;
sol = ode45(odefunc, tspan, z0, odeOptions);
timeConstraint = toc;

% Solve ODE with stiff spring (problem 20 method) for comparison
odefuncSpring = @(t, z) twoMassesSpringRHS(t, z, p);
tic;
solSpring = ode45(odefuncSpring, tspan, z0, odeOptions);
timeSpring = toc;

disp("Constraint ode45 time = " + num2str(timeConstraint) + " s, steps = " + num2str(length(sol.x)));
disp("Spring (k = " + num2str(p.k) + ") ode45 time = " + num2str(timeSpring) + " s, steps = " + num2str(length(solSpring.x)));

% Unpack ODE sol
tarray = linspace(t0, tf, 1e3+1);
zarray = deval(sol, tarray);

x1array = zarray(1, :);
v1array = zarray(2, :);
x2array = zarray(3, :);
v2array = zarray(4, :);
larray = x2array - x1array;


%%% Arrays for F_G vs F0
% should your numerical solution always have the  
% property that F = m_tot a_G same as F0
xcomarray = (m1 * x1array + m2 * x2array) / (m1 + m2);
vcomarray = gradient(xcomarray, tarray);
acomarray = gradient(vcomarray, tarray);
F_G = (m1 + m2) * acomarray;

diffInF = F_G - F0;
diffInF = diffInF(3:end-2);
tarrayForDiffInF = tarray(3:end-2);

rmsErrorInFGF0 = sqrt(mean(diffInF .^ 2));
rmsErrorInFGF0Array = rmsErrorInFGF0 * ones(size(tarrayForDiffInF));

%%% Arrays for satisfaction of constraint
% constraint = x2 - x1 = d
errorInConstraintArray = (x2array - x1array) - d;
rmsErrorInConstraint = sqrt(mean(errorInConstraintArray.^2));

%%% Arrays for constraint vs spring comparison
zarraySpring = deval(solSpring, tarray);
x1arraySpring = zarraySpring(1, :);
x2arraySpring = zarraySpring(3, :);

absErrorX1 = abs(x1arraySpring - x1array);
absErrorX2 = abs(x2arraySpring - x2array);
absErrorInConstraintSpring = abs((x2arraySpring - x1arraySpring) - d);

rmsErrorX1 = sqrt(mean(absErrorX1.^2));
rmsErrorX2 = sqrt(mean(absErrorX2.^2));
rmsErrorInConstraintSpring = sqrt(mean(absErrorInConstraintSpring.^2));

%%% Arrays for COM accuracy vs tolerance
% exact COM: (m1 + m2) xG'' = Fext, starting at rest
xcomExact = xcomarray(1) + integral(@(s) (tf - s) * p.F0_func(s), t0, tf, 'ArrayValued', true) / (m1 + m2);

tolArray = logspace(-3, -12, 10);
comErrorConstraintArray = zeros(size(tolArray));
comErrorSpringArray = zeros(size(tolArray));

for i = 1:length(tolArray)
    odeOptionsTol = odeset('RelTol', tolArray(i), 'AbsTol', tolArray(i));

    solTol = ode45(odefunc, tspan, z0, odeOptionsTol);
    zendTol = deval(solTol, tf);
    xcomEndTol = (m1 * zendTol(1) + m2 * zendTol(3)) / (m1 + m2);
    comErrorConstraintArray(i) = abs(xcomEndTol - xcomExact);

    solSpringTol = ode45(odefuncSpring, tspan, z0, odeOptionsTol);
    zendSpringTol = deval(solSpringTol, tf);
    xcomEndSpringTol = (m1 * zendSpringTol(1) + m2 * zendSpringTol(3)) / (m1 + m2);
    comErrorSpringArray(i) = abs(xcomEndSpringTol - xcomExact);
end

% Plot
figure('Name', 'Trajectory')
tiledlayout(1, 1);

nexttile;
plot(tarray, x1array, 'r-'); hold on;
plot(tarray, x2array, 'b-');
plot(tarray, xcomarray, 'k--');

legend('x_1', 'x_2', 'x_{COM}');
xlabel('t (s)');
ylabel('position (m)');
grid on;

figure('Name', 'ErrorInConstraint')
plot(tarray, errorInConstraintArray, 'DisplayName', 'Error'); hold on;
plot(tarray, rmsErrorInConstraint * ones(size(tarray)), '-.', 'DisplayName', "RMS Error = " + num2str(rmsErrorInConstraint)); 
title('Error in constrain $= (x_2 - x_1) - d$', 'Interpreter', 'latex')
grid on;
legend();

figure('Name', 'ConstraintVsSpring')
tiledlayout(3, 1);

nexttile;
semilogy(tarray, absErrorX1, 'r-', 'DisplayName', "RMS = " + num2str(rmsErrorX1));
title('$|x_1^{spring} - x_1^{constraint}|$', 'Interpreter', 'latex')
grid on;
legend();

nexttile;
semilogy(tarray, absErrorX2, 'b-', 'DisplayName', "RMS = " + num2str(rmsErrorX2));
title('$|x_2^{spring} - x_2^{constraint}|$', 'Interpreter', 'latex')
grid on;
legend();

nexttile;
semilogy(tarray, absErrorInConstraintSpring, 'k-', 'DisplayName', "Spring, RMS = " + num2str(rmsErrorInConstraintSpring)); hold on;
semilogy(tarray, abs(errorInConstraintArray), 'g-', 'DisplayName', "Constraint, RMS = " + num2str(rmsErrorInConstraint));
title('$|(x_2 - x_1) - d|$', 'Interpreter', 'latex')
xlabel('t (s)');
grid on;
legend();

figure('Name', 'COMErrorVsTolerance')
loglog(tolArray, comErrorConstraintArray, 'g.-', 'DisplayName', 'Constraint'); hold on;
loglog(tolArray, comErrorSpringArray, 'k.-', 'DisplayName', "Spring, k = " + num2str(p.k));
title('Error in $x_G(t_f)$ vs ode45 tolerance', 'Interpreter', 'latex')
xlabel('RelTol = AbsTol');
ylabel('|x_G(t_f) - x_G^{exact}(t_f)|');
grid on;
legend();


% Animation
figure('Name', 'Animation');

rodLine = plot([x1array(1) x2array(1)], [0 0], 'k-', 'LineWidth', 2); hold on;
mass1Data = plot(x1array(1), 0, 'r.', 'MarkerSize', 40 * sqrt(m1));
mass2Data = plot(x2array(1), 0, 'b.', 'MarkerSize', 40 * sqrt(m2));
comData = plot(xcomarray(1), 0, 'kx', 'MarkerSize', 12);

legend('Rod', 'm1', 'm2', 'COM');
xlabel('x (m)');
grid on;
axis([min([x1array x2array]) - 0.1, max([x1array x2array]) + 0.1, -1, 1]);
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
        rodLine.XData = [x1now x2now];
        mass1Data.XData = x1now;
        mass2Data.XData = x2now;
        comData.XData = xcomnow;
    end

    drawnow limitrate;
    t = toc * speedupFactor;
end
