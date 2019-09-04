%% getGMD
% Calculates the GMD of a coil to itself given the width and height of the coil.
% All distances in meters.

function GMD = getGMD(h, w)
    d = sqrt(h^2 + w^2);
    u = w^2/h^2*log(d^2/w^2);
    v = h^2/w^2*log(d^2/h^2);
    x = w/h*atan(h/w);
    y = h/2*atan(w/h);
    theta = (u + v + 25)/12 - 2/3*(x + y);
    GMD = d*exp(-theta);
end