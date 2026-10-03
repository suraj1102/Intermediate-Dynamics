
function zdot = cannonBallRHS(t, z, p)
    % Anchor / Fixed hinge at origin
    [m, g, c] = unpack_p(p);
    
    i_hat = [1; 0]; j_hat = [0; 1];

    r_vec = z(1:2);
    v_vec = z(3:4);

    rdot_vec = v_vec;

    F_weight = m*g* -j_hat;
    F_damping = c * -v_vec;

    F_tot = F_damping + F_weight;
    vdot_vec = F_tot ./ m;

    zdot = [rdot_vec; vdot_vec];
end
