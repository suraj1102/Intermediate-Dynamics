function zdot = myRHS(t, z, p)
    x = z(1:2);
    v = z(3:4);

    i_hat = [1 0];
    j_hat = [0 1];
    
    F1 = sin(t) * i_hat + cos(t) * j_hat;
    F2 = cos(t) * i_hat + sin(t) * j_hat;

    xdot = v;
    vdot = F1 + F2;

    zdot = [xdot vdot];
end
