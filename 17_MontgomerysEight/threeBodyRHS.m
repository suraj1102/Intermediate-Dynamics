function zdot = threeBodyRHS(t, z, p)
    r1 = z(1:2);
    r2 = z(3:4);
    r3 = z(5:6);
    v1 = z(7:8);
    v2 = z(9:10);
    v3 = z(11:12);

    G = p.G; m1 = p.m1; m2 = p.m2; m3 = p.m3;

    r12 = r2 - r1;
    r13 = r3 - r1;
    r23 = r3 - r2;

    F12 = G * m1 * m2 / norm(r12)^2 * r12 ./ norm(r12);
    F13 = G * m1 * m3 / norm(r13)^2 * r13 ./ norm(r13);
    F23 = G * m2 * m3 / norm(r23)^2 * r23 ./ norm(r23);

    v1dot = (F12 + F13) ./ m1;
    v2dot = (-F12 + F23) ./ m2;
    v3dot = (-F13 - F23) ./ m3;

    zdot = [v1; v2; v3; v1dot; v2dot; v3dot];
end
