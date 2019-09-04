%% computePERMS
% Calculates the RMS value of a waveform using the power electronics convention
% of separating DC and AC components and computing the square root of the sum of
% their squares.

function rmsVal = computePERMS(wf)
    DC = mean(wf);
    AC = wf - DC; % variation about DC
    rmsVal = sqrt(DC^2 + rms(AC)^2);
end