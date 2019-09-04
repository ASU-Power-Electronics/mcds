%% m2AWG
% AWG of a wire from cross-sectional area in meters squared.

function AWG = m2AWG(A)
    AWG = 36 - 19.5.*log(A.*1e6./0.012668)./log(92);
end