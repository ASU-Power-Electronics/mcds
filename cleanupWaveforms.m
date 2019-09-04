%% cleanupWaveforms
% Checks for and removes zero-length time steps prior to interpolation.

function [w, t] = cleanupWaveforms(wf, tv, T)
    dtVec = diff([tv; tv(1) + T]);
    idx = [];
    
    for i = 1:length(dtVec) - 1
        if dtVec(i) > eps
            idx = [idx; i];
        end
    end
    idx = idx + ones(length(idx), 1);
    
    t = [tv(1); tv(idx); T];
    w = [wf(1, :); wf(idx, :); wf(1, :)];
end