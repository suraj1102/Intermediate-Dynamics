function drawCircle(origin, radius)
    xorigin = origin(1);
    yorigin = origin(2);

    theta = linspace(0, 2*pi, 1e2+1);

    xarray = xorigin + radius .* cos(theta);
    yarray = yorigin + radius .* sin(theta);

    plot(xarray, yarray, 'k--');
end