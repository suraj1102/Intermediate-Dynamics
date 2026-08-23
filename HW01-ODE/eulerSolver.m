function [t, z] = eulerSolver(odefunc, tspan, z0, h)
    t = tspan(1) : h : tspan(2);
    t = t'; % column vector
    z = zeros(size(t, 1), length(z0)); % N, state_dim

    z(1, :) = z0;

    for i = 1 : size(z, 1) - 1
        t_now = t(i);
        z_now = z(i, :);
        z_next = z_now + h * odefunc(t_now, z_now);

        z(i+1, :) = z_next;
    end
end
