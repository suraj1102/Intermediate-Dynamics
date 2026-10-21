function zdot = doublePendDAErhs(t, z, p)
    % Z = [rE vE rP vP]

    m1 = p.m1; m2 = p.m2;
    l1 = p.l1; l2 = p.l2;
    g = p.g;

    rE = z(1:2); vE = z(3:4);
    rP = z(5:6); vP = z(7:8);

    x1 = rE(1); y1 = rE(2);
    x2 = rP(1); y2 = rP(2);
    x1dot = vE(1); y1dot = vE(2);
    x2dot = vP(1); y2dot = vP(2);

    rPE = rP - rE;
    
    % Unit Vectors
    iHat = [1 0]';
    jHat = [0 1]';
    rEHat = rE / norm(rE);
    rPEHat = rPE / norm(rPE);

    elbowConstraint = -x1dot^2 - y1dot^2;
    forearmConstraint = -(  x1dot^2 + x2dot^2 + y1dot^2 + y2dot^2 + ...
        -2*x1dot*x2dot + -2*y1dot*y2dot);
    F = [0 -m1*g 0 -m2*g elbowConstraint forearmConstraint]';

    A = [
        m1 0 0 0 -dot(rEHat, iHat) -dot(rPEHat, iHat);
        0 m1 0 0 -dot(rEHat, jHat) -dot(rPEHat, jHat);
        0 0 m2 0  0                 dot(rPEHat, iHat);
        0 0  0 m2 0                 dot(rPEHat, jHat);
        x1 y1 0 0 0 0;
        (x1-x2) (y1-y2) (x2-x1) (y2-y1) 0 0
    ];

    u = A\F;
    aE = u(1:2);
    aP = u(3:4);
    zdot = [vE; aE; vP; aP];
end