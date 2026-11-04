clc; close all; clear;

% theta
syms t
theta = (1 - cos(t))*pi;
dtheta = diff(theta, t);
ddtheta = diff(dtheta, t);

% variables
syms phi(t) m
dphi = diff(phi, t);
ddphi = diff(dphi, t);

% unit vectors
EX = [1; 0; 0];
EY = [0; 1; 0];
k  = [0; 0; 1];

ex = cos(phi) * EX + sin(phi) * EY;
ey = -sin(phi) * EX + cos(phi) * EY;

erB = cos(theta) * ex + sin(theta) * ey;
etheta = -sin(theta) * ex + cos(theta) * ey;

% omega
omega = dphi * k;
domega = ddphi * k;

% position vectors
r1_O = 2 * ey;
r2_O = -2 * ey;

rB = 2 * ex;
radiusOfRotatingMotion = 1;
r3_O = rB + radiusOfRotatingMotion * erB;

% AMB of mass 1
a1_O = cross(domega, r1_O) + cross(omega, cross(omega, r1_O));
dH1_O = cross(r1_O, m * a1_O);

a1_0_direct = diff(r1_O, t, 2);
dH1_O_direct = cross(r1_O, m * a1_0_direct);


syms p dp ddp
eq1 = subs(dH1_O, [diff(phi,t,2), diff(phi,t), phi], [ddp, dp, p]);
simplify(eq1)

eq2 = subs(dH1_O_direct, [diff(phi,t,2), diff(phi,t), phi], [ddp, dp, p]);
simplify(eq2)
