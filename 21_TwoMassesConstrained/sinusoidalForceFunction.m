function Fext = sinusoidalForceFunction(t, p)
    F0 = p.F0;
    Fext = F0 * sin(t);
end