
function zdot = springMassRHS(t, z, p)
    % Anchor / Fixed hinge at origin
    [m, g, k, l0] = unpack_p(p);
    
    i_hat = [1; 0]; j_hat = [0; 1];

    r_vec = z(1:2);
    v_vec = z(3:4);

    rdot_vec = v_vec;


    F_weight = m*g* -j_hat;
    
    l = norm(r_vec);
    lambda_hat = r_vec / l;
    F_s = -k * (l - l0) * lambda_hat;

    F_tot = F_s + F_weight;
    vdot_vec = F_tot ./ m;

    zdot = [rdot_vec; vdot_vec];
end
