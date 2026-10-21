function zdot = twoMassesSpringRHS(t, z, p)
    x1 = z(1); v1 = z(2);
    x2 = z(3); v2 = z(4);
    m1 = p.m1; m2 = p.m2;
    k  = p.k;  l0 = p.l0;
    F0 = p.F0;

    l = (x2-x1);
    dl = (l-l0);
    Fs = dl * k;

    a1 = Fs/m1;
    a2 = (-Fs + F0)/m2;

    zdot = [v1; a1; v2; a2];
end