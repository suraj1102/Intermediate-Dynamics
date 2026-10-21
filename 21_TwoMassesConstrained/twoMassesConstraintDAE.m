function zdot = twoMassesSpringRHS(t, z, p)
    x1 = z(1); v1 = z(2);
    x2 = z(3); v2 = z(4);
    m1 = p.m1; m2 = p.m2;
    F0_func = p.F0_func;

    Fext = F0_func(t);

    A = [m1 0 -1;
        0 m2 1;
        -1 1 0];
    F = [0 Fext 0]';
    u = A\F;

    a = u(1:2);
    a1 = a(1);
    a2 = a(2);
    
    zdot = [v1; a1; v2; a2];
end