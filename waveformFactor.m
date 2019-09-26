%% waveformFactor
% Voltage waveform factor calculation.  `v` is primary voltage, and symmetric is
% a boolean value representing whether or not B is symmetric wrt to zero.

function Kf = waveformFactor(v, t, symmetric)
    % compute flux to find zero-to-peak time
    phi = zeros(size(v));
    dt = t(2) - t(1);
    
    for time = 2:length(phi)
        phi(time) = phi(time - 1) + v(time)*dt;
    end
    
    if symmetric
        phi = phi - mean(phi);
    end
    
%     figure
%     plot(t, phi)
%     grid on
    
    % find max
    [phiMax, idx2] = max(phi(:));
    
    % find first zero-crossing and rising cycle after it
    idx1 = find(abs(phi) <= phiMax/10, 1);
    [~, idx1] = min(abs(phi(1:idx1 + 256)));
    
    % set zero-to-peak time region
    tau = t(idx1:idx2);
    
    % create voltage mask
    trash = zeros(1, idx1 - 1);
    comp = [trash, tau];
    comp = [comp, zeros(1, length(t) - length(comp))];
    
    if symmetric
        comp(1) = inf; % always finds first element otherwise
    end
    
    mask = t == comp;
    
    % compute denominator values
    tau0 = tau(end) - tau(1);
    T0 = t(end) - t(1);
    
    % get zero-to-peak flux region in voltage waveform
    v_k = v.*mask;
    
%     figure
%     plot(t, v_k)
%     grid on
    
    % compute form factor of waveform (k = V_RMS/V_avg(in tau))
    VRMS = rms(v);
    Vavg = mean(v_k(idx1:idx2));
    k = abs(VRMS/Vavg);
    
    % compute waveform factor
    Kf = k*(T0/tau0);
end