%% AWG2m
% Diameter and cross-sectional area of an AWG wire in meters and meters squared.
% Accepts both scalar and vector inputs and returns a vector or matrix, as
% appropriate.

function [d, A] = AWG2m(awg)
    d = 0.127.*92.^((36 - awg)./39).*1e-3;
    A = pi/4.*d.^2;
end