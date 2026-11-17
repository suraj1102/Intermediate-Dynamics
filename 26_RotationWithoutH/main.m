%% main.m  -  Rotation with zero angular momentum (Problem 26)
% Solves the EoM for the plate angle phi(t) with ode45, plots the results
% and the trajectories of the masses, and animates the plate as seen from
% the fixed frame F.
%
% Notation (same as the write-up):
%   F      fixed frame with axes X', Y' (origin O')
%   O      plate frame, unit vectors ex, ey; O == O' for all time
%   phi    rotation of the plate relative to F (CCW positive)
%   theta  angle of m3 on its circle about B, measured from the plate x axis
%   m1, m2 welded at y = +2, -2;  m3 on circle about B = (2,0);  m4 = -m3

clc; close all; clear;

%% Parameters
d          = 2;      % |OB| = |OA| = distance of m1, m2 from O
R          = 1;      % radius of the circles
Lplate     = 4.2;    % half-width of the drawn plate (visual only)
nCycles    = 3;      % number of revolutions of m3/m4 (try 3 to see the ratchet)
nPerCycle  = 300;    % output points (= animation frames) per revolution
saveGif    = false;  % true -> also write the animation to a GIF
gifName    = 'zero_H_rotation.gif';

%% EoM from the AMB:  phi'' = f(t, phi')   (valid for 0 <= t <= pi)
ddphi = @(t, dp) -(pi*(2*cos(t).*cos(pi*cos(t)) - cos(t) + 2*pi*sin(pi*cos(t)) ...
        + 4*dp.*sin(t).*sin(pi*cos(t)) - 2*pi*cos(t).^2.*sin(pi*cos(t)))) ...
        ./ (4*cos(pi*cos(t)) - 9);

%% Solve with ode45, one revolution at a time
% State z = [phi; phi'].  Each revolution uses the local time
% tau = t - (c-1)*pi in [0, pi], so theta goes 0 -> 2*pi every cycle.
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
t = [];  Z = [];
z = [0; 0];                                   % starts from rest
for c = 1:nCycles
    t0 = (c - 1)*pi;
    tc = linspace(t0, t0 + pi, nPerCycle).';
    [tk, Zk] = ode45(@(tt, zz) [zz(2); ddphi(tt - t0, zz(2))], tc, z, opts);
    z = Zk(end, :).';
    if c > 1                                  % drop duplicated boundary point
        tk(1) = [];  Zk(1, :) = [];
    end
    t = [t; tk];  Z = [Z; Zk];                %#ok<AGROW>
end
phi  = Z(:, 1);
dphi = Z(:, 2);
N    = numel(t);

% prescribed motion (unwrapped over cycles)
cyc  = min(floor(t/pi), nCycles - 1);
tau  = t - cyc*pi;
th   = pi*(1 - cos(tau)) + 2*pi*cyc;
dth  = pi*sin(tau);

%% Results
phiExact = nCycles*pi*(7/sqrt(65) - 1);
fprintf('Net rotation after %d revolution(s):\n', nCycles);
fprintf('  ode45       : phi = %+.8f rad  (%+.4f deg)\n', phi(end), rad2deg(phi(end)));
fprintf('  closed form : phi = %+.8f rad  (%+.4f deg)\n', phiExact, rad2deg(phiExact));

%% Positions and velocities in F (2 x N arrays, columns = time)
ex = [ cos(phi).';  sin(phi).'];
ey = [-sin(phi).';  cos(phi).'];
cth = cos(th).';  sth = sin(th).';

r1 =  d*ey;
r2 = -d*ey;
r3 = (d + R*cth).*ex + (R*sth).*ey;
r4 = -r3;

% velocities: d/dt ex = phi' ey,  d/dt ey = -phi' ex
v1 = -d*dphi.'.*ex;
v2 = -v1;
v3 = (-R*sth.*(dth.' + dphi.')).*ex + (R*cth.*dth.' + (d + R*cth).*dphi.').*ey;
v4 = -v3;

% angular momentum about O per unit mass (should stay 0)
crossz = @(r, v) r(1, :).*v(2, :) - r(2, :).*v(1, :);
Hz = crossz(r1, v1) + crossz(r2, v2) + crossz(r3, v3) + crossz(r4, v4);
fprintf('  max |H_z/m| over the motion = %.2e\n', max(abs(Hz)));

%% Plots: time histories
figure('Name', 'Time histories', 'Color', 'w', 'Position', [80 80 1000 650]);

subplot(2, 2, 1);
plot(t, phi, 'LineWidth', 1.5); hold on; grid on;
yline(phiExact, '--r', 'closed form');
xlabel('t'); ylabel('\phi [rad]'); title('Plate rotation \phi(t)');

subplot(2, 2, 2);
plot(t, dphi, 'LineWidth', 1.5); grid on;
xlabel('t'); ylabel('d\phi/dt [rad/s]'); title('Plate angular velocity');

subplot(2, 2, 3);
plot(t, th, 'LineWidth', 1.5); grid on;
xlabel('t'); ylabel('\theta [rad]'); title('Prescribed \theta(t) of m_3, m_4');

subplot(2, 2, 4);
plot(t, Hz, 'LineWidth', 1.5); grid on;
xlabel('t'); ylabel('H_z / m'); title('Angular momentum check (should be 0)');

%% Plots: trajectories
figure('Name', 'Trajectories', 'Color', 'w', 'Position', [120 120 1100 520]);

% (a) seen from the fixed frame F
subplot(1, 2, 1); hold on; axis equal; grid on; box on;
plot(r1(1, :), r1(2, :), 'k', 'LineWidth', 1.5);
plot(r2(1, :), r2(2, :), 'k', 'LineWidth', 1.5);
plot(r3(1, :), r3(2, :), 'r', 'LineWidth', 1.5);
plot(r4(1, :), r4(2, :), 'b', 'LineWidth', 1.5);
Pstart = [r1(:, 1) r2(:, 1) r3(:, 1) r4(:, 1)];
Pend   = [r1(:, N) r2(:, N) r3(:, N) r4(:, N)];
plot(Pstart(1, :), Pstart(2, :), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
plot(Pend(1, :),   Pend(2, :),   'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
plot([-Lplate Lplate], [0 0], 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
plot([0 0], [-Lplate Lplate], 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
plot(Lplate*[-ex(1, N) ex(1, N)], Lplate*[-ex(2, N) ex(2, N)], 'k:');
plot(Lplate*[-ey(1, N) ey(1, N)], Lplate*[-ey(2, N) ey(2, N)], 'k:');
xlabel('X'''); ylabel('Y''');
title(sprintf('Paths in F  (final plate axes dotted, \\phi = %.1f%s)', ...
      rad2deg(phi(end)), char(176)));
legend({'m_1', 'm_2', 'm_3', 'm_4'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
axis([-1 1 -1 1]*(Lplate + 0.5));

% (b) seen from the plate frame O
subplot(1, 2, 2); hold on; axis equal; grid on; box on;
q = linspace(0, 2*pi, 200);
plot(0, d, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(0, -d, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(d + R*cos(q),  R*sin(q), 'r', 'LineWidth', 1.5);
plot(-d - R*cos(q), -R*sin(q), 'b', 'LineWidth', 1.5);
plot([-Lplate Lplate], [0 0], 'k--');
plot([0 0], [-Lplate Lplate], 'k--');
xlabel('x'); ylabel('y');
title('Same paths seen from the plate');
axis([-1 1 -1 1]*(Lplate + 0.5));

%% Animation (seen from F)
figA = figure('Name', 'Animation', 'Color', 'w', 'Position', [160 60 700 720]);
axA = axes('Parent', figA); hold(axA, 'on'); axis(axA, 'equal');
lim = Lplate*sqrt(2) + 0.3;
axis(axA, [-lim lim -lim lim]); box(axA, 'on');
xlabel(axA, 'X'''); ylabel(axA, 'Y''');

rot = @(a) [cos(a) -sin(a); sin(a) cos(a)];

% fixed axes of F
plot(axA, [-lim lim], [0 0], '--', 'Color', [0.7 0.7 0.7]);
plot(axA, [0 0], [-lim lim], '--', 'Color', [0.7 0.7 0.7]);
text(axA, lim - 0.5, 0.25, 'X''', 'Color', [0.5 0.5 0.5]);
text(axA, 0.15, lim - 0.3, 'Y''', 'Color', [0.5 0.5 0.5]);

% geometry in plate coordinates
plateP  = Lplate*[-1 1 1 -1; -1 -1 1 1];
xAxisP  = [-Lplate Lplate; 0 0];
yAxisP  = [0 0; -Lplate Lplate];
circB_P = [d + R*cos(q);  R*sin(q)];
circA_P = [-d + R*cos(q); R*sin(q)];

% graphics objects
hPlate = patch(axA, plateP(1, :), plateP(2, :), [0.93 0.93 0.93], ...
               'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.5);
hXax   = plot(axA, xAxisP(1, :), xAxisP(2, :), 'k-', 'LineWidth', 1);
hYax   = plot(axA, yAxisP(1, :), yAxisP(2, :), 'k-', 'LineWidth', 1);
hCB    = plot(axA, circB_P(1, :), circB_P(2, :), ':', 'Color', [0.8 0.3 0.3]);
hCA    = plot(axA, circA_P(1, :), circA_P(2, :), ':', 'Color', [0.3 0.3 0.8]);
hLine  = plot(axA, [0 0], [0 0], 'k:', 'LineWidth', 1);            % m3-O-m4
hTr3   = plot(axA, NaN, NaN, 'r-', 'LineWidth', 1);                % trails
hTr4   = plot(axA, NaN, NaN, 'b-', 'LineWidth', 1);
hTr1   = plot(axA, NaN, NaN, 'k-', 'LineWidth', 1);
hTr2   = plot(axA, NaN, NaN, 'k-', 'LineWidth', 1);
hM12   = plot(axA, [0 0], [0 0], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
hM3    = plot(axA, 0, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
hM4    = plot(axA, 0, 0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 12);
plot(axA, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
hTxt   = text(axA, -lim + 0.3, lim - 0.4, '', 'FontSize', 11, 'FontName', 'Courier');
title(axA, 'Free plate, H = 0: masses m_3, m_4 circle, plate turns');

for kk = 1:N
    Q  = rot(phi(kk));
    Pp = Q*plateP;   Xa = Q*xAxisP;   Ya = Q*yAxisP;
    CB = Q*circB_P;  CA = Q*circA_P;

    set(hPlate, 'XData', Pp(1, :), 'YData', Pp(2, :));
    set(hXax,   'XData', Xa(1, :), 'YData', Xa(2, :));
    set(hYax,   'XData', Ya(1, :), 'YData', Ya(2, :));
    set(hCB,    'XData', CB(1, :), 'YData', CB(2, :));
    set(hCA,    'XData', CA(1, :), 'YData', CA(2, :));
    set(hLine,  'XData', [r3(1, kk) r4(1, kk)], 'YData', [r3(2, kk) r4(2, kk)]);
    set(hTr3,   'XData', r3(1, 1:kk), 'YData', r3(2, 1:kk));
    set(hTr4,   'XData', r4(1, 1:kk), 'YData', r4(2, 1:kk));
    set(hTr1,   'XData', r1(1, 1:kk), 'YData', r1(2, 1:kk));
    set(hTr2,   'XData', r2(1, 1:kk), 'YData', r2(2, 1:kk));
    set(hM12,   'XData', [r1(1, kk) r2(1, kk)], 'YData', [r1(2, kk) r2(2, kk)]);
    set(hM3,    'XData', r3(1, kk), 'YData', r3(2, kk));
    set(hM4,    'XData', r4(1, kk), 'YData', r4(2, kk));
    set(hTxt, 'String', sprintf('t     = %5.2f\ntheta = %6.1f deg\nphi   = %6.2f deg', ...
        t(kk), rad2deg(th(kk)), rad2deg(phi(kk))));
    drawnow;

    if saveGif
        frame = getframe(figA);
        [imind, cmap] = rgb2ind(frame2im(frame), 256);
        if kk == 1
            imwrite(imind, cmap, gifName, 'gif', 'LoopCount', inf, 'DelayTime', 0.03);
        else
            imwrite(imind, cmap, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
        end
    else
        pause(0.01);
    end
end
