function zdot = inverseGravityRHS(t, z, p);
    r1 = z(1:2);
    r2 = z(3:4);
    v1 = z(5:6);
    v2 = z(7:8);

    G = p.G; m1 = p.m1; m2 = p.m2;

    r = r2 - r1;
    r_hat = r ./ norm(r);

    F = G * m1 * m2 / norm(r)^2 * r_hat;

    v1dot = F ./ m1;
    v2dot = -F ./ m2;

    zdot = [v1; v2; v1dot; v2dot];
end