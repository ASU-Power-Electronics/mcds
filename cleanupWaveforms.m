%% cleanupWaveforms
% Checks for and perturbs near-zero-length time steps prior to interpolation.

function [w, t] = cleanupWaveforms(wf, tv, T)
    dtVec = diff([tv; tv(1) + T]);
%     idx = [];
    
    idxZeroPre = dtVec <= eps;
    idxZeroPost = circshift(idxZeroPre, 1);
    tv(idxZeroPre) = tv(idxZeroPre) - T/1e6;
    tv(idxZeroPost) = tv(idxZeroPost) + T/1e6;
    
    if ~isequal(tv(1), 0)
        tv(1) = 0;
    end
    
    t = [tv(1:end - 1); T];
    
    w = wf;
    
%     for tPos = 1:length(dtVec) - 1
%         if dtVec(tPos) > eps
%             idx = [idx; tPos];
%         end
%     end
%     
%     idx = idx + ones(length(idx), 1);
%     
%     t = [tv(1); tv(idx); T];
%     w = [wf(1, :); wf(idx, :); wf(1, :)];
end