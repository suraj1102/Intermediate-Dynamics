clc; close all; clear;

syms t m
syms phi(t)            % phi is a FUNCTION of t, so diff(phi,t) is not 0
ph = phi(t);           % plain sym expression (avoids symfun indexing pitfalls)

% prescribed motion of the x-axis masses
theta = (1 - cos(t))*pi;

% fixed unit vectors
EX = [1; 0; 0];  EY = [0; 1; 0];  k = [0; 0; 1];

% plate unit vectors
ex = cos(ph)*EX + sin(ph)*EY;
ey = -sin(ph)*EX + cos(ph)*EY;
erB = cos(theta)*ex + sin(theta)*ey;

% positions relative to O (O is fixed, O = O')
r1 = 2*ey;
r2 = -2*ey;
r3 = 2*ex + erB;
r4 = -r3;
R  = {r1, r2, r3, r4};

% inertial accelerations: just differentiate twice in the fixed frame
dH = sym(zeros(3,1));
for i = 1:4
    a_i = diff(R{i}, t, 2);
    dH  = dH + cross(R{i}, m*a_i);
end

% check: mass 1 alone should give [0; 0; 4*m*phi'']
dH1 = simplify(cross(r1, m*diff(r1, t, 2)))

% AMB: dH_z = 0 -> solve for phi''
syms p dp ddp
eq = subs(dH(3), [diff(ph,t,2), diff(ph,t), ph], [ddp, dp, p]);
ddphi_expr = simplify(solve(eq == 0, ddp))

% integrate phi'' numerically, starting from rest
f = matlabFunction(ddphi_expr, 'Vars', {t, p, dp});
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[tt, Y] = ode45(@(tt, y) [y(2); f(tt, y(1), y(2))], [0 pi], [0; 0], opts);

fprintf('phi(pi) = %.6f rad = %.2f deg\n', Y(end,1), rad2deg(Y(end,1)));
fprintf('exact   = %.6f rad\n', pi*(7/sqrt(65) - 1));

plot(tt, Y(:,1), 'LineWidth', 1.5); grid on
xlabel('t'); ylabel('\phi(t) [rad]')