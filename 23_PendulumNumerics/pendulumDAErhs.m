function zdot = pendulumDAErhs(t, z, p)
    m = p.m; g = p.g;
    r = z(1:2);
    v = z(3:4);
    l = norm(r);

    A = [
        m 0 r(1)/l;
        0 m r(2)/l;
        r(1) r(2) 0
    ];
    F = [0 -m*g -v(1)^2-v(2)^2]';

    u = A\F;
    a = u(1:2);

    zdot = [v; a];
end