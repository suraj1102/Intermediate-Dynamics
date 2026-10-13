function zdot = myRHS(t, z, p)
    zdot = zeros(length(z), 1);
    G = p.G; m = [p.m1 p.m2 p.m3];

    for i = 1:2:5
        x_i = z(i: i+1);
        v_i = z(i+6: i+6+1);
        m_i = m((i + 1) / 2);

        ii = mod(i + 1, 6) + 1;
        iii = mod(i + 3, 6) + 1;

        x_ii = z(ii: ii+1);
        x_iii = z(iii: iii+1);

        d1 = -(x_i - x_ii);
        d2 = -(x_i - x_iii);

        F1 = G * m_i * m((ii + 1) / 2) * d1 / norm(d1)^3;
        F2 = G * m_i * m((iii + 1) / 2) * d2 / norm(d2)^3;

        Ftot = F1 + F2;
        a_i = Ftot ./ m_i;

        zdot(i:i+1) = v_i;
        zdot(i+6:i+6+1) = a_i;
    end
end