%% calculateFlux
% Pointwise flux integrator.  Input v is the turns-scaled voltage waveform.
% Implementation of v = Ndphi/dt due to Faraday's Law.

function phi = calculateFlux(v, t)
%     dt = t(2) - t(1);
%     g = @(t) v(round(t/dt) + ones(1, length(t)));
%     phi = zeros(size(v));
    
%     for i = 1:length(v) - 1
%         if v(i) % v ~= 0 -> non-zero slope
%             phi(i) = integral(g, t(1), t(i), 'Waypoints', t(1:i));
%         else % avoid rounding/truncation error
%             if i == 1
%                 phi(i) = 0; % dc value assumed not known, symmetry shift later
%             else            
%                 phi(i) = phi(i - 1);
%             end
%         end
%     end

    phi = cumtrapz(t, v);
    
    phi = phi - mean(phi);
    phi(end) = phi(1); % maintain periodicity
end