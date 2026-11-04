close all; clear; clc;

p = struct();
p.numPoints = 1e3;

% z = [phi, dphi]
z0 = [0, 0];
tspan = [0 pi];

rhsFunc = @(t, z) myRHS(t, z, p);
sol = ode45(rhsFunc, tspan, z0);


tArray = linspace(tspan(1), tspan(2), p.numPoints+1);
zArray = deval(sol, tArray);

phiArray = zArray(1, :);
dphiArray = zArray(2, :);

% Plot

figure()
tiledlayout(2, 1);

nexttile;
plot(tArray, phiArray)
xlabel("t"); ylabel("\phi")
grid on;

nexttile;
plot(tArray, dphiArray)
xlabel("t"); ylabel("$\dot{\phi}$", Interpreter="latex")
grid on;







function zdot = myRHS(t, z, p)
    theta_t = @(t) (1-cos(t)) * pi;
    dtheta_t = @(t)pi.*sin(t);
    ddtheta_t = @(t)pi.*cos(t);

    theta = theta_t(t);
    dtheta = dtheta_t(t);
    ddtheta = ddtheta_t(t);


    % Unpack z = [phi, dphi]

    phi = z(1);
    dphi = z(2);



    t1 = 8 * dtheta * dphi;
    t2 = -2 * (2*theta - 1) * dphi^2;
    t3 = -4 * theta * dtheta * dphi;
    t4 = 2 * theta * dtheta;
    t5 = -4 * dtheta;
    t6 = -6 * ddtheta;

    ddphi = (t1+t2+t3+t4+t5+t6) / (26 - 4*theta);


    zdot = [];
    zdot(1) = dphi;
    zdot(2) = ddphi;
    zdot = zdot';
end