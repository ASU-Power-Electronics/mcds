%% airInductance
% Calculates mutual inductance of two coils (or a coil with itself) given the
% geometry of the coil(s), neglecting contribution of any core.  Returns two
% values, base and leakage inductance.

function [Lab, Lal] = airInductance(N1, N2, A, R, length, GMD, tol)
    global MU_0
    
    sum = 0;
    k = 1;
    delta = 1;
    
    % infinite image sum with truncation
    while delta > tol
        betaK = 2*pi*k/length;
        % next 2 terms listed in reverse in literature (A<->R)
        b1 = besseli(1, betaK*A);
        b2 = besselk(1, betaK*R);
        c1 = cos(betaK*GMD);
        old = sum;
        sum = sum + b1*b2*c1; % accumulated leakage flux
        delta = abs(old - sum);
        k = k + 1;
        if isnan(sum) || isinf(sum)
%             fprintf('No good:  air inductance outta whack at k = %d.\n', k)
            sum = old; % put it back
            break
        end
    end
    
    Lab = MU_0*N1*N2*A*R*2*pi/length*R/(2*A);
    Lal = MU_0*N1*N2*A*R*2*pi/length*2*sum;
end