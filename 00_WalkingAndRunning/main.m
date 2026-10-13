clear; close all; clc;

%% Settings
gait = 'trajOpt';   % 'myICs', 'flatWalk', 'run' or 'trajOpt'

% System Params
systemParams.m = 1;
systemParams.g = 1;
systemParams.lmax = 1; % max leg extension - optim constraint

m = systemParams.m;  g = systemParams.g;

% System Variables
systemVariables.xc = 0;
systemVariables.tStanceStart = 0;
systemVariables.numPoints = 1e3+1;
systemVariables.ode45Options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

switch gait
    case 'myICs'
        systemParams.ts = 1;
        systemParams.tf = 0;
        systemParams.d  = 1;
        systemVariables.F = @(tInStance) constantLegForce(tInStance);

        x0 = 0;
        y0 = 0.6;
        vx0 = 1;
        vy0 = 0.05;

    case 'flatWalk'
        % Body stays at height h during stance; no flight phase
        h = 0.8;
        systemParams.ts = 1;
        systemParams.tf = 0;
        systemParams.d  = 1;
        w = sqrt(g/h);
        A = (systemParams.d/2) / sinh(w*systemParams.ts/2);
        ts_ = systemParams.ts;
        systemVariables.F = @(tau) flatWalkForce(tau, m, g, h, A, w, ts_);

        x0 = -0.5;  y0 = 0.8;  vx0 = A*w*cosh(w*ts_/2);  vy0 = 0;   % vx0 = 1.102059

    case 'run'
        % Sine-hump leg force, symmetric about mid-stance
        Fpk = 3;
        systemParams.ts = 0.5;
        systemParams.tf = 0.438697;
        systemParams.d  = 1.458185;
        ts_ = systemParams.ts;
        systemVariables.F = @(tau) sineLegForce(tau, Fpk, ts_);

        x0 = -0.383918;  y0 = 0.893941;  vx0 = 1.573636;  vy0 = -0.219348;

    case 'trajOpt'
        systemParams.v_des = 1;
        systemParams.d  = 1.0;
        systemParams.tStep = systemParams.d / systemParams.v_des;
        systemParams.collocationPoints = 10;
        systemParams.Fmax = 100;

        [Fgrid, x0, y0, vx0, vy0, tStance] = performTrajOpt(systemParams, systemVariables);

        systemParams.ts = tStance;
        systemParams.tf = systemParams.tStep - tStance;
        tGrid = linspace(0, tStance, numel(Fgrid));
        systemVariables.F = @(tau) interp1(tGrid, Fgrid, min(max(tau, 0), tStance));
end

z0 = [x0; y0; vx0; vy0];
systemVariables.xc = z0(1) - x0;   % = 0: foot at origin for step 1



%% Integrate - Multi Step
numSteps = 10;
tArray = [];
zArray = [];

current_z0 = z0;
tStart = 0;
for step = 1:numSteps
    [tArrayStep, zArrayStep] = oneStep(current_z0, tStart, systemParams, systemVariables);

    % From step 2 on, the first sample repeats the previous step's last sample
    if step > 1
        tArrayStep = tArrayStep(2:end);
        zArrayStep = zArrayStep(:, 2:end);
    end

    tArray = [tArray tArrayStep];
    zArray = [zArray zArrayStep];

    % Update things from one step to next
    current_z0 = zArrayStep(:, end);
    systemVariables.xc = current_z0(1) - x0;
    tStart = tArray(end);
end


%% Plot
ts = systemParams.ts;  tf = systemParams.tf;
tStep = ts + tf;
N = numel(tArray);
assert(size(zArray, 2) == N, 'tArray and zArray have different lengths');

% Recover step number and phase from time (valid because ts, tf are fixed)
tol = 1e-9;
stepArray = min(floor((tArray + tol) / tStep) + 1, numSteps);
tInStep   = tArray - (stepArray - 1) * tStep;
isStance  = tInStep <= ts + tol;

% Foot for each step = body x at the start of that step (your current foot rule)
xcArray = zeros(1, numSteps);
for s = 1:numSteps
    xcArray(s) = zArray(1, find(stepArray == s, 1)) - x0;
end

cStance = [0 0.4470 0.7410];       % blue
cFlight = [0.8500 0.3250 0.0980];  % orange

figure('Name', 'Multi step', 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%d steps:  t_s = %.2f,  t_f = %.2f', numSteps, ts, tf));

axPath = nexttile(tl, [1 2]); hold(axPath, 'on');
axVel  = nexttile(tl, [1 2]); hold(axVel,  'on');
axT = gobjects(1, 4);
for k = 1:4
    axT(k) = nexttile(tl); hold(axT(k), 'on');
end

% Draw each (step, phase) segment; extend by one sample so segments join
for s = 1:numSteps
    for phase = [true false]
        idx = find(stepArray == s & isStance == phase);
        if isempty(idx), continue; end
        if idx(end) < N, idx = [idx, idx(end) + 1]; end

        if phase, c = cStance; else, c = cFlight; end
        plot(axPath, zArray(1, idx), zArray(2, idx), 'Color', c, 'LineWidth', 1.5);
        plot(axVel,  zArray(3, idx), zArray(4, idx), 'Color', c, 'LineWidth', 1.5);
        for k = 1:4
            plot(axT(k), tArray(idx), zArray(k, idx), 'Color', c, 'LineWidth', 1.5);
        end
    end
end

% Body path
plot(axPath, xcArray, zeros(1, numSteps), 'k^', 'MarkerFaceColor', 'k');  % feet
plot(axPath, zArray(1,1), zArray(2,1), 'ko', 'MarkerFaceColor', 'w');     % start
yline(axPath, 0, 'k-');                                                   % ground
axis(axPath, 'equal'); grid(axPath, 'on');
xlabel(axPath, 'x'); ylabel(axPath, 'y'); title(axPath, 'Body path');
hS = plot(axPath, NaN, NaN, 'Color', cStance, 'LineWidth', 1.5);
hF = plot(axPath, NaN, NaN, 'Color', cFlight, 'LineWidth', 1.5);
hP = plot(axPath, NaN, NaN, 'k^', 'MarkerFaceColor', 'k');
legend(axPath, [hS hF hP], {'Stance', 'Flight', 'Feet'}, 'Location', 'best');

% Velocity
plot(axVel, zArray(3,1), zArray(4,1), 'ko', 'MarkerFaceColor', 'w');
grid(axVel, 'on');
xlabel(axVel, 'v_x'); ylabel(axVel, 'v_y'); title(axVel, 'Velocity');

% States vs time
tLiftoff   = (0:numSteps-1) * tStep + ts;   % stance -> flight
tTouchdown = (1:numSteps-1) * tStep;        % flight -> next stance
labels = {'x', 'y', 'v_x', 'v_y'};
for k = 1:4
    for tl_ = tLiftoff,   xline(axT(k), tl_, 'k--'); end
    for td_ = tTouchdown, xline(axT(k), td_, 'k:');  end
    grid(axT(k), 'on');
    xlabel(axT(k), 't'); ylabel(axT(k), labels{k}); title(axT(k), [labels{k} '(t)']);
end

%% Animate
playbackSpeed = 1;      % simulation time units per real second
viewHalfWidth = 1.5;    % half-width of the window that follows the body
tfinal = tArray(end);
lmax = systemParams.lmax;

figAnim = figure('Name', 'Animation', 'Color', 'w');
axA = axes(figAnim); hold(axA, 'on'); grid(axA, 'on');
daspect(axA, [1 1 1]);
xlabel(axA, 'x'); ylabel(axA, 'y');
ylim(axA, [-0.1 1.2]);

% Static elements
plot(axA, zArray(1,:), zArray(2,:), ':', 'Color', [0.75 0.75 0.75]);        % full path, faint
plot(axA, xcArray, zeros(size(xcArray)), 'k^', 'MarkerFaceColor', 'k');    % footholds
yline(axA, 0, 'k-', 'LineWidth', 1.5);                                     % ground

% Dynamic elements
hTrail = plot(axA, NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
hLeg   = plot(axA, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 3);
hBody  = plot(axA, NaN, NaN, 'o', 'MarkerSize', 14, ...
              'MarkerFaceColor', cStance, 'MarkerEdgeColor', 'k');
hText  = text(axA, 0.02, 0.96, '', 'Units', 'normalized', ...
              'VerticalAlignment', 'top', 'FontName', 'FixedWidth');

% Animation loop
i = 1;          % index into tArray; only moves forward
tic;
tSim = 0;
while tSim < tfinal
    if ~isgraphics(figAnim), break; end      % stop cleanly if the window is closed

    tSim = toc * playbackSpeed;

    % Fetch nearest point: advance i to the last sample <= tSim, then pick the closer neighbour
    while i < N && tArray(i+1) <= tSim
        i = i + 1;
    end
    j = i;
    if i < N && abs(tArray(i+1) - tSim) < abs(tArray(i) - tSim)
        j = i + 1;
    end

    x = zArray(1, j);  y = zArray(2, j);
    s = stepArray(j);

    % Update plots
    if isStance(j)
        xc = xcArray(s);
        l  = hypot(x - xc, y);
        set(hLeg, 'XData', [xc x], 'YData', [0 y], 'Visible', 'on');
        set(hBody, 'MarkerFaceColor', cStance);
        phaseStr = sprintf('stance   l = %.3f / %.2f', l, lmax);
    else
        set(hLeg, 'Visible', 'off');
        set(hBody, 'MarkerFaceColor', cFlight);
        phaseStr = 'flight';
    end

    set(hBody,  'XData', x, 'YData', y);
    set(hTrail, 'XData', zArray(1, 1:j), 'YData', zArray(2, 1:j));
    xlim(axA, x + [-1 1] * viewHalfWidth);
    hText.String = sprintf('t = %5.2f   step %d   %s', tArray(j), s, phaseStr);

    drawnow limitrate
end


%% Function Definitions
function [Fgrid, x0, y0, vx0, vy0, tStance] = performTrajOpt(systemParams, systemVariables)
    % Decision vector z = [x0, y0, vx0, vy0, F_0 ... F_N, ts]  (length N+6)
    N = systemParams.collocationPoints;
    m = systemParams.m;  g = systemParams.g;
    lmax = systemParams.lmax;  Fmax = systemParams.Fmax;
    d = systemParams.d;  tStep = systemParams.tStep;

    % Initial guess: flat walk (no flight) sampled on the grid
    h = 0.8;
    tsGuess = tStep;
    w = sqrt(g/h);
    A = (d/2) / sinh(w*tsGuess/2);
    tauGrid = linspace(0, tsGuess, N+1);
    FGuess = arrayfun(@(tau) flatWalkForce(tau, m, g, h, A, w, tsGuess), tauGrid);
    guessArray = [-d/2, h, A*w*cosh(w*tsGuess/2), 0, FGuess, tsGuess];

    % Sanity check on the guess before optimizing
    CGuess = trajOptCost(guessArray, systemParams, systemVariables);
    [cGuess, ceqGuess] = trajOptConstraints(guessArray, systemParams, systemVariables);
    fprintf('Guess: C = %.5f (flat walk d/(8h) = %.5f)\n', CGuess, d/(8*h));
    fprintf('Guess: |ceq| = %.2e,  max(c) = %.3f\n', norm(ceqGuess), max(cGuess));

    % Bounds
    lb = [-lmax, 0,    -Inf, -Inf, zeros(1, N+1),      0.01];
    ub = [ lmax, lmax,  Inf,  Inf, Fmax*ones(1, N+1),  tStep];

    nonlcon = @(z) trajOptConstraints(z, systemParams, systemVariables);
    cost = @(z) trajOptCost(z, systemParams, systemVariables);
    % oneStep integrates the whole stance with one adaptive ode45 call, and the kinks in the
    % piecewise-linear F make its output slightly noisy. The default finite-difference step
    % (~1e-8) turns that noise into wrong gradients, so use a bigger, central step.
    options = optimoptions('fmincon', 'Algorithm','sqp', ...
        'Display','iter', 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e3, ...
        'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-5);

    [outputArray, fval, exitflag, output] = fmincon(cost, guessArray, ...
        [], [], [], [], lb, ub, nonlcon, options);

    [~, ceqOpt] = trajOptConstraints(outputArray, systemParams, systemVariables);
    fprintf('fmincon exitflag %d after %d iterations: C = %.5f, |ceq| = %.2e\n', ...
        exitflag, output.iterations, fval, norm(ceqOpt));

    x0 = outputArray(1);  y0 = outputArray(2);
    vx0 = outputArray(3); vy0 = outputArray(4);
    Fgrid = outputArray(5:N+5);
    tStance = outputArray(N+6);

    fprintf('Optimum: ts = %.4f, tf = %.4f\n', tStance, tStep - tStance);
    fprintf('Optimum: x0 = %.4f, y0 = %.4f, vx0 = %.4f, vy0 = %.4f\n', x0, y0, vx0, vy0);
    fprintf('Optimum: F grid = [%s]\n', sprintf(' %.3f', Fgrid));
end

function [c, ceq] = trajOptConstraints(z, systemParams, systemVariables)
    [z0, systemParams, systemVariables, tGrid] = unpackTrajOpt(z, systemParams, systemVariables);
    [tArray, zArray] = oneStep(z0, 0, systemParams, systemVariables);

    % Periodicity: same state one stride later
    zEnd = zArray(:, end);
    ceq = [zEnd(1) - z0(1) - systemParams.d;
           zEnd(2) - z0(2);
           zEnd(3) - z0(3);
           zEnd(4) - z0(4)];

    % Leg can't be longer than lmax at any grid point
    % (oneStep returns numPoints stance samples first, then flight)
    iStance = 1:systemVariables.numPoints;
    lStance = hypot(zArray(1, iStance) - systemVariables.xc, zArray(2, iStance));
    lGrid = interp1(tArray(iStance), lStance, tGrid);
    c = lGrid(:) - systemParams.lmax;
end

function cost = trajOptCost(z, systemParams, systemVariables)
    [z0, systemParams, systemVariables] = unpackTrajOpt(z, systemParams, systemVariables);
    [tArray, zArray] = oneStep(z0, 0, systemParams, systemVariables);

    % Leg power during stance (leg is massless, so it does no work in flight)
    iStance = 1:systemVariables.numPoints;
    t = tArray(iStance);
    x = zArray(1, iStance) - systemVariables.xc;  y = zArray(2, iStance);
    vx = zArray(3, iStance);                      vy = zArray(4, iStance);
    l = sqrt(x.^2 + y.^2);
    ldot = (x.*vx + y.*vy) ./ l;
    F = arrayfun(systemVariables.F, t);
    power = F .* ldot;

    % Only positive work costs energy; negative work is free
    work = trapz(t, max(power, 0));

    % Cost of transport: positive work / (m g d)
    cost = work / (systemParams.m * systemParams.g * systemParams.d);
end

function [z0, systemParams, systemVariables, tGrid] = unpackTrajOpt(z, systemParams, systemVariables)
    % Turn the decision vector into the inputs oneStep expects
    N = systemParams.collocationPoints;
    Fgrid = z(5:N+5);
    tStance = z(N+6);
    tGrid = linspace(0, tStance, N+1);

    systemParams.ts = tStance;
    systemParams.tf = systemParams.tStep - tStance;
    systemVariables.F = @(tau) interp1(tGrid, Fgrid, min(max(tau, 0), tStance));
    systemVariables.xc = 0;   % foot at origin, so z(1) = x0 is the body x relative to the foot

    z0 = [z(1); z(2); z(3); z(4)];
end



function [tArray, zArray] = oneStep(z0, tStart, systemParams, systemVariables)
    ts = systemParams.ts; tf = systemParams.tf;
    numPoints = systemVariables.numPoints;
    systemVariables.tStanceStart = tStart;

    % Stance
    tSpanStance = [tStart tStart+ts];
    stanceODEHandle = @(t, z) stanceODE(t, z, systemParams, systemVariables);
    stanceSol = ode45(stanceODEHandle, tSpanStance, z0, systemVariables.ode45Options);

    tArray = linspace(tSpanStance(1), tSpanStance(2), numPoints);
    zArray = deval(stanceSol, tArray);

    % Flight (skipped for gaits with no flight phase)
    if tf > 0
        stateAtEndOfStance = zArray(:, end);
        tSpanFlight = [tStart+ts tStart+ts+tf];
        flightODEHandle = @(t, z) flightODE(t, z, systemParams, systemVariables);
        flightSol = ode45(flightODEHandle, tSpanFlight, stateAtEndOfStance, systemVariables.ode45Options);

        tFlightArray = linspace(tSpanFlight(1), tSpanFlight(2), numPoints);
        tFlightArray = tFlightArray(2:end); % Delete duplicate point
        zArray = [zArray deval(flightSol, tFlightArray)];
        tArray = [tArray tFlightArray];
    end
end

function F = constantLegForce(tInStance)
    F = 1.5;
end


function zdot = stanceODE(t, z, systemParams, systemVariables)
    x = z(1); y = z(2);
    vx = z(3); vy = z(4);
    
    m = systemParams.m;
    g = systemParams.g;

    xc = systemVariables.xc;
    F = systemVariables.F(t - systemVariables.tStanceStart);

    l = sqrt((x - xc)^2 + y^2);

    vxdot = F * (x - xc) / (l * m);
    vydot = (-m*g + F*y/l) / m;

    zdot = [vx; vy; vxdot; vydot];
end


function zdot = flightODE(t, z, systemParams, systemVariables)
    x = z(1); y = z(2);
    vx = z(3); vy = z(4);
    
    m = systemParams.m;
    g = systemParams.g;

    vxdot = 0;
    vydot = -g;

    zdot = [vx; vy; vxdot; vydot];
end


function F = flatWalkForce(tau, m, g, h, A, w, ts)
    % Keeps the body at constant height h: F*y/l = m*g with x(tau) = A*sinh(w*(tau - ts/2))
    x = A * sinh(w * (tau - ts/2));
    F = m * g * sqrt(x^2 + h^2) / h;
end

function F = sineLegForce(tau, Fpk, ts)
    F = Fpk * sin(pi * tau / ts);
end
