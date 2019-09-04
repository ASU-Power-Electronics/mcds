%% coreInductance
% Calculates mutual inductance of two coils (or a coil with itself) given the
% geometry of the coil(s), due to the contribution of the core.  Returns two
% values, base and leakage inductance.

function [Lcb, Lcl] = coreInductance(N1, N2, rIn1, rOut1, rIn2, rOut2, h1, w1, h2, w2, length, b, mu_i, sigma, omega)
    global MU_0
    
    % necessary functions (uses Struve01 from Mathworks)
    p = @(z) pi*z*(besselk(1, z)*StruveL0(z) + StruveL1(z)*besselk(0, z))/2;
    P = @(x, y, beta) (p(beta*x) - p(beta*y))/beta^2;
    Q = @(x, y, beta) 2*(cos(beta*(x - y)/2) - cos(beta*(x + y)/2))/beta^2;
    G = @(beta) sqrt(beta^2 + 1j*omega*mu_i*MU_0*sigma); % gamma
    fNum = @(beta) 1 - (besseli(1, beta*b)*G(beta)*b*besseli(0, G(beta)*b))/(mu_i*beta*b*besseli(0, beta*b)*besseli(1, G(beta)*b));
    fDen = @(beta) 1 + (besselk(1, beta*b)*G(beta)*b*besseli(0, G(beta)*b))/(mu_i*beta*b*besselk(0, beta*b)*besseli(1, G(beta)*b));
    f = @(beta) fNum(beta)/fDen(beta);
    phi = @(beta) besseli(0, beta*b)/besselk(0, beta*b)*real(f(beta));
    
    sum = 0;
    k = 1;
    delta = 1;
    
    % infinite image sum with truncation
    while delta > 1e3*eps
        betaK = 2*pi*k/length;
        P1 = P(rOut1, rIn1, betaK);
        P2 = P(rOut2, rIn2, betaK);
        QVal = Q(w1, w2, betaK);
        phiVal = phi(betaK);
        old = sum;
        sum = sum + P1*P2*QVal*phiVal;
        delta = abs(old - sum);
        k = k + 1;
        if isnan(sum) || isinf(sum)
            fprintf('No good:  inductance outta whack.\n')
            sum = old; % put it back
            break
        end
    end
    
    Lcb = mu_i*MU_0*N1*N2*pi*b^2/length;
    Lcl = MU_0*N1*N2*2*pi*2*sum/(h1*w1*h2*w2*length);
end