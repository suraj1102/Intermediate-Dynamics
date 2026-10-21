function zdot = pendulumLMBrhs(t, z, p)
    % FUNCTION z = [theta thetadot] NOT [x y vx vy]
    
    l = p.l; g = p.g;

    theta = z(1);
    thetadot = z(2);

    thetaddot = -g/l * cos(theta);

    zdot = [thetadot; thetaddot];
end