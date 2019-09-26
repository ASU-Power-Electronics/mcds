%% Niner
% Creates a reduced-length vector capturing large spikes in the waveform
% provided.  w0 and wp0 are the waveform and its derivative in time, and t0 is
% the associated time vector.  The name comes from the original test case, from
% which nine time points were extracted to sufficiently model the current given
% a trapezoidal waveform with dead time.

function [t, idx] = niner(t0, w0, wp0)    
    [~, sortIdx] = sort(diff(wp0));
    idx = zeros(1, 2);
    nPoints = 4;
    
    % for a quality 9-point signal, difference between points should be close to
    % length/8 on average; if not, reduce points to spread them out
    while abs(128 - mean(diff(idx))) > 32 && nPoints > 0
        idx = sort([sortIdx(1:nPoints), sortIdx(end - (nPoints - 1):end)]);
        nPoints = nPoints - 1;
    end
    
    % take actual data instead of difference location
    idx = idx + ones(size(idx));
    
    
    % include first and last point in reduced vector to ensure periodicity
    if ~(idx(1) == 1)
        idx = [1, idx];
    end
    
    if ~(idx(end) == length(t0))
        idx = [idx, length(t0)];
    end
    
    t = t0(idx);
    w = w0(idx);
%     w0RMS = rms(w0(1:end - 1));
%     
%     figure
%     hold on
%     plot(t0, w0)
%     plot(t, w)
%     hold off
%     grid on
%     refline(0, w0RMS)
%     title('Ideal Winding Current Waveform & Linear Approximation')
%     xlabel('t (s)')
%     ylabel('i_w (A)')
%     legend('Ideal', 'PWL Approximation', 'RMS Value')
end