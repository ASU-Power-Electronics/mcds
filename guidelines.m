%% Guidelines
% Computes litz wire sizing guidelines.
%TODO: Add option to suggest parallel wires for large current densities.

function result = guidelines(P, S, b)
    global J_MAX
    
    [~, nwp] = size(P);
    [~, nws] = size(S);
    result = struct();
    result.P = struct(); % results for primary winding(s)
    result.S = struct(); % results for secondary winding(s)
    
    % Litz wire strand gauges and areas
    awgz = 32:48;
    [~, Az] = AWG2m(awgz);
    
    % available wire constructions (eBay) (note, no 34AWG wire for sale)
    constructions = {[19, 50]; ...
                     0; ...
                     [16, 27, 41, 65, 105, 165, 265]; ...
                     [5, 10, 15, 16, 20, 25, 30, 32, 40, 45, 50, 60, 64, 66, 70, 80, 90, 100, 128, 140, 160, 180, 200, 220, 250, 280, 320, 350, 400, 420, 450, 540, 600, 640, 660, 700, 800, 900, 1000, 1050, 1300, 1650, 2000, 2500]; ...
                     [17, 27, 42, 66, 100, 105, 108, 170, 270, 435, 700]; ...
                     [16, 26, 40, 66, 105, 165, 270, 420, 660, 1050]; ...
                     [10, 20, 30, 40, 60, 100, 160, 170, 255, 270, 405, 420, 650, 660, 1050]; ...
                     [10, 20, 30, 40, 60, 66, 100, 175, 180, 220, 250, 270, 330, 420, 660, 924, 1162, 1200]; ...
                     [180, 300, 420, 675, 1100]};
    
    % table from Sullivan (2014)
	k = [130, 203, 318, 496, 771, 1.2e3, 1.8e3, 2.8e3, 4.4e3, 6.7e3, 10e3, 16e3, 24e3, 36e3, 54e3, 79e3, 115e3];
	d = [0.202, 0.180, 0.160, 0.143, 0.127, 0.113, 0.101, 0.090, 0.080, 0.071, 0.063, 0.056, 0.050, 0.045, 0.040, 0.035, 0.032];
    
    if nwp > 1
        for p = 1:nwp
            N = P(p).N;
            delta = P(p).delta_p*1e3;
            
            result.P(p).Amin = P(p).I_pRMS/J_MAX; % resulting minimum cross-sectional area of conductor
            result.P(p).maxAWG = floor(m2AWG(P(p).A_pMax)); % theoretical maximum strand size, AWG
            result.P(p).minbAWG = ceil(m2AWG(result.P(p).Amin)); % minimum bundle size, AWG
            [result.P(p).d, ~] = AWG2m(result.P(p).minbAWG); % minimum bundle diameter
            result.P(p).nS = k.*delta^2*b/N; % number of strands for AWG sizes 32-48
            result.P(p).D = sqrt(result.P(p).nS).*d*1e-3/0.68125; % associated outer diameter
            result.P(p).A = Az.*result.P(p).nS; % associated area
            junk = P(p).d_pMax - d.*1e-3; % placeholder for difference
            [~, idx] = find(junk > 0, 1); % index of closest safe fit
            result.P(p).dMax = d(idx)*1e-3; % maximum strand diameter (skin effect)
            result.P(p).minAWG = awgz(idx); % resulting strand size, AWG
            result.P(p).ndMax = ceil(k(idx).*delta^2*b/N); % associated number of strands
            result.P(p).Db = sqrt(result.P(p).ndMax).*result.P(p).dMax/0.68125; % bundle outer diameters (corrected for litz)
            result.P(p).nbpl = floor(b*1e-3/result.P(p).Db); % associated bundles per layer (to fit in bobbin breadth)
            cons = constructions{idx} - result.P(p).ndMax; % purchasable constructions for AWG
            idx2 = find(cons >= 0, 1); % identify all number of strands greater than or equal to
            result.P(p).sugStr = [num2str(cons(idx2) + result.P(p).ndMax), '/', num2str(awgz(idx))]; % formatted string for suggestion
        end
    else
        N = P.N;
        delta = P.delta_p*1e3;

        result.P.Amin = P.I_pRMS/J_MAX;
        result.P.maxAWG = floor(m2AWG(P.A_pMax));
        result.P.minbAWG = ceil(m2AWG(result.P.Amin));
        [result.P.d, ~] = AWG2m(result.P.minbAWG);
        result.P.nS = k.*delta^2*b/N;
        result.P.D = sqrt(result.P.nS).*d*1e-3/0.68125;
        result.P.A = Az.*result.P.nS;
        junk = P.d_pMax - d.*1e-3;
        [~, idx] = find(junk > 0, 1);
        result.P.dMax = d(idx)*1e-3;
        result.P.minAWG = awgz(idx);
        result.P.ndMax = ceil(k(idx).*delta^2*b/N);
        result.P.Db = sqrt(result.P.ndMax).*result.P.dMax/0.68125;
        result.P.nbpl = floor(b*1e-3/result.P.Db);
        cons = constructions{idx} - result.P.ndMax;
        idx2 = find(cons >= 0, 1);
        result.P.sugStr = [num2str(cons(idx2) + result.P.ndMax), '/', num2str(awgz(idx))];
    end
    
    if nws > 1
        for s = 1:nws
            N = S(s).N;
            delta = S(s).delta_s*1e3;
            
            result.S(s).Amin = S(s).I_sRMS/J_MAX;
            result.S(s).maxAWG = floor(m2AWG(S(s).A_sMax));
            result.S(s).minbAWG = ceil(m2AWG(result.S(s).Amin));
            [result.S(s).d, ~] = AWG2m(result.S(s).minbAWG);
            result.S(s).nS = k.*delta^2*b/N;
            result.S(s).D = sqrt(result.S(s).nS).*d*1e-3/0.68125;
            result.S(s).A = Az.*result.S(s).nS;
            junk = S(s).d_sMax - d.*1e-3;
            [~, idx] = find(junk > 0, 1);
            result.S(s).dMax = d(idx)*1e-3;
            result.S(s).minAWG = awgz(idx);
            result.S(s).ndMax = ceil(k(idx).*delta^2*b/N);
            result.S(s).Db = sqrt(result.S(s).ndMax).*result.S(s).dMax/0.68125;
            result.S(s).nbpl = floor(b*1e-3/result.S(s).Db);
            cons = constructions{idx} - result.S(s).ndMax;
            idx2 = find(cons >= 0, 1);
            result.S(s).sugStr = [num2str(cons(idx2) + result.S(s).ndMax), '/', num2str(awgz(idx))];
        end
    else
        N = S.N;
        delta = S.delta_s*1e3;

        result.S.Amin = S.I_sRMS/J_MAX;
        result.S.maxAWG = floor(m2AWG(S.A_sMax));
        result.S.minbAWG = ceil(m2AWG(result.S.Amin));
        [result.S.d, ~] = AWG2m(result.S.minbAWG);
        result.S.nS = k.*delta^2*b/N;
        result.S.D = sqrt(result.S.nS).*d*1e-3/0.68125;
        junk = S.d_sMax - d.*1e-3;
        [~, idx] = find(junk > 0, 1);
        result.S.dMax = d(idx)*1e-3;
        result.S.minAWG = awgz(idx);
        result.S.ndMax = ceil(k(idx).*delta^2*b/N);
        result.S.Db = sqrt(result.S.ndMax).*result.S.dMax/0.68125;
        result.S.nbpl = floor(b*1e-3/result.S.Db);
        cons = constructions{idx} - result.S.ndMax;
        idx2 = find(cons >= 0, 1);
        result.S.sugStr = [num2str(cons(idx2) + result.S.ndMax), '/', num2str(awgz(idx))];
    end
end