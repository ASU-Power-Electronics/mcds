%% evaluatePermeability
% Computes relative permeability given a material and peak flux density value.
% If material.mu_a is a scalar, that value is returned.  Otherwise, the value is
% computed by an n-term Gaussian sum, where n is one-third the number of
% elements in mu_a.

function mu_r = evaluatePermeability(mu_a, B_pk)
    if isequal(numel(mu_a), 1)
        mu_r = mu_a;
    else
        mu_r = 0;
        
        for idx = 1:3:numel(mu_a)
            an = mu_a(idx);
            bn = mu_a(idx + 1);
            cn = mu_a(idx +2);
            mu_r = mu_r + an*exp(-((B_pk - bn)/cn)^2);
        end
    end
end