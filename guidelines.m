%% Guidelines
% Computes litz wire sizing guidelines.
%TODO: Add option to suggest parallel wires for large current densities.

function result = guidelines(P, S, b, h)
    global J_MAX
    
    [~, nwp] = size(P);
    [~, nws] = size(S);
    result = struct();
    result.P = struct(); % results for primary winding(s)
    result.P.constructions = struct();
    result.S = struct(); % results for secondary winding(s)
    result.S.constructions = struct();
    result.Wgap = 0;
    
    resp = inputdlg('Height of additional space/insulation between windings [mm]:', ...
                    'Adjust inter-winding space', ...
                    1, {'0'});
	result.Wgap = str2double(resp)*1e-3;
    h = h - (nwp + nws - 1)*result.Wgap;
    
    Ab = h*b; % bobbin window area [m]^2
    ACutot = 0;
    
    % reduced proximity loss, high surface area to volume ratio Sullivan (2014)
    Jtarget = J_MAX*(nwp + nws);
    
    % Litz wire strand gauges and areas
    awgz = 32:2:48;
    [~, Az] = AWG2m(awgz);
    
    % available wire constructions (eBay) (note, no 34AWG wire for sale)
    load('Wires.mat', 'Wires')
    
    % table from Sullivan paper (2014) - k is a parameter, d is strand diameter
    % truncated by removing odd gauges, since they aren't purchasable
	k = [130, 318, 771, 1.8e3, 4.4e3, 10e3, 24e3, 54e3, 115e3]*1e9;
	d = [0.202, 0.160, 0.127, 0.101, 0.080, 0.063, 0.050, 0.040, 0.032]*1e-3;
    
    % find minimum areas for Jtarget and total minimum area of all windings
    if nwp > 1
        for p = 1:nwp
            result.P(p).Amin = max(abs(P(p).waveform.i_p))/Jtarget;
            ACutot = ACutot + result.P(p).Amin*P(p).N;
        end
    else
        result.P.Amin = max(abs(P.waveform.i_p))/Jtarget;
        ACutot = ACutot + result.P.Amin*P.N;
    end
    
    if nws > 1
        for s = 1:nws
            result.S(s).Amin = max(abs(S(s).waveform.i_s))/Jtarget;
            ACutot = ACutot + result.S(s).Amin*S(s).N;
        end
    else
        result.S.Amin = max(abs(S.waveform.i_s))/Jtarget;
        ACutot = ACutot + result.S.Amin*S.N;
    end
    
    if nwp > 1
        for p = 1:nwp
            thisW = P(p);
            thisR = result.P(p);
            
            % structure arrays need fields pre-defined
            result.P(p).Amax = [];
            result.P(p).AWGmin = [];
            result.P(p).AWGmax = [];
            result.P(p) = computeResultsP(thisW, thisR);
        end
    else
        thisW = P;
        thisR = result.P;
        result.P = computeResultsP(thisW, thisR);
    end
    
    if nws > 1
        for s = 1:nws
            thisW = S(s);
            thisR = result.S(s);
            
            % structure arrays need fields pre-defined
            result.S(s).Amax = [];
            result.S(s).AWGmin = [];
            result.S(s).AWGmax = [];
            result.S(s) = computeResultsS(thisW, thisR);
        end
    else
        thisW = S;
        thisR = result.S;
        result.S = computeResultsS(thisW, thisR);
    end
    
    %% Nested Functions
    
    % computeResultsP
    % Computes results for primary.    
    function res = computeResultsP(W, res)
        validCons = 0;       % valid construction counter
        N = W.N;
        delta = W.delta_p;
        
        % Skin depth provides range of bundles
        dmax = W.d_pMax;
        [~, idx] = find(d(dmax > d), 1);    % index of closest safe fit
        ArangeS = Az(idx:end);              % strand copper areas
        AWGrangeS = awgz(idx:end);          % AWG values
        
        % base case
        nS = round(k(idx:end)*delta^2*b/N); % integer number of strands
        ArangeB = ArangeS.*nS;              % bundle copper areas for nS
        AWGrangeB = ceil(m2AWG(ArangeB));	% equivalent integer bundle AWGs
        [DeffB, ~] = AWG2m(AWGrangeB);      % bundle effective Cu diameters
        
        % +25% number of strands
        nShi = ceil(nS*1.25);
        ArangeBhi = ArangeS.*nShi;
        AWGrangeBhi = ceil(m2AWG(ArangeBhi));
        [DeffBhi, ~] = AWG2m(AWGrangeBhi);
        
        % -25% number of strands
        nSlo = floor(nS*0.75);
        ArangeBlo = ArangeS.*nSlo;
        AWGrangeBlo = ceil(m2AWG(ArangeBlo));
        [DeffBlo, ~] = AWG2m(AWGrangeBlo);
        
        % Winding area bounds bundle size above
        AmaxTot = Ab*N*res.Amin/ACutot;	% fraction of bobbin window
        AmaxOut = (pi/4)*AmaxTot/N;     % max outer round bundle area for fit
        
        % Final upper bound is minimum of window and skin depth max
        res.Amax = AmaxOut*(0.68125)^2; % max bundle copper area
        res.AWGmax = m2AWG(res.Amax);   % max bundle effective AWG
        
        % Current density bounds bundle size below
        res.AWGmin = m2AWG(res.Amin);   % min bundle effective AWG
        
        figure
        plot(ArangeB)
        hold on
        plot(ArangeBhi)
        plot(ArangeBlo)
        refline(0, res.Amin)
        refline(0, res.Amax)
        hold off
        
        for awg = idx:length(AWGrangeS) + idx - 1
            isValid = Wires.AWG(awg).A_b >= res.Amin & ...
                      Wires.AWG(awg).A_b <= res.Amax & ...
                      Wires.AWG(awg).N_s > nSlo(awg) & ...
                      Wires.AWG(awg).N_s < nShi(awg);
            
            if any(isValid)
                for validCon = 1:length(isValid)                    
                    if isValid(validCon)
                        validCons = validCons + 1;
                        construction.N_s = Wires.AWG(awg).N_s(validCon);
                        construction.AWG = Wires.AWGvals(awg);
                        construction.AWG_e = round(Wires.AWG(awg).AWG_e(validCon));
                        construction.A_b = Wires.AWG(awg).A_b(validCon);
                        construction.D_o = Wires.AWG(awg).D_o(validCon);
                        construction.N_bpl = floor(b./Wires.AWG(awg).D_o(validCon));
                        
                        if isempty(fieldnames(res.constructions))
                            res.constructions.N_s = [];
                            res.constructions.AWG = [];
                            res.constructions.AWG_e = [];
                            res.constructions.A_b = [];
                            res.constructions.D_o = [];
                            res.constructions.N_bpl = [];
                        end
                        
                        res.constructions(validCons) = construction;
                    end
                end
            end
            
            clear isValid
        end
    end
    
    % computeResultsS
    % Computes results for secondary.
    function res = computeResultsS(W, res)
        validCons = 0;       % valid construction counter
        N = W.N;
        delta = W.delta_s;
        
        % Skin depth provides range of bundles
        dmax = W.d_sMax;
        [~, idx] = find(d(dmax > d), 1);    % index of closest safe fit
        ArangeS = Az(idx:end);              % strand copper areas
        AWGrangeS = awgz(idx:end);          % AWG values
        
        % base case
        nS = round(k(idx:end)*delta^2*b/N); % integer number of strands
        ArangeB = ArangeS.*nS;              % bundle copper areas for nS
        AWGrangeB = ceil(m2AWG(ArangeB));	% equivalent integer bundle AWGs
        [DeffB, ~] = AWG2m(AWGrangeB);      % bundle effective Cu diameters
        
        % +25% number of strands
        nShi = ceil(nS*1.25);
        ArangeBhi = ArangeS.*nShi;
        AWGrangeBhi = ceil(m2AWG(ArangeBhi));
        [DeffBhi, ~] = AWG2m(AWGrangeBhi);
        
        % -25% number of strands
        nSlo = floor(nS*0.75);
        ArangeBlo = ArangeS.*nSlo;
        AWGrangeBlo = ceil(m2AWG(ArangeBlo));
        [DeffBlo, ~] = AWG2m(AWGrangeBlo);
        
        % Winding area bounds bundle size above
        AmaxTot = Ab*N*res.Amin/ACutot;	% fraction of bobbin window
        AmaxOut = (pi/4)*AmaxTot/N;     % max outer round bundle area for fit
        res.Amax = AmaxOut*(0.68125)^2; % max bundle copper area
        res.AWGmax = m2AWG(res.Amax);   % max bundle effective AWG
        
        % Current density bounds bundle size below
        res.AWGmin = m2AWG(res.Amin);   % min bundle effective AWG
        
        figure
        plot(ArangeB)
        hold on
        plot(ArangeBhi)
        plot(ArangeBlo)
        refline(0, res.Amin)
        refline(0, res.Amax)
        hold off
        
        for awg = idx:length(AWGrangeS) + idx - 1
            isValid = Wires.AWG(awg).A_b >= res.Amin & ...
                      Wires.AWG(awg).A_b <= res.Amax & ...
                      Wires.AWG(awg).N_s > nSlo(awg) & ...
                      Wires.AWG(awg).N_s < nShi(awg);
            
            if any(isValid)
                for validCon = 1:length(isValid)                    
                    if isValid(validCon)
                        validCons = validCons + 1;
                        construction.N_s = Wires.AWG(awg).N_s(validCon);
                        construction.AWG = Wires.AWGvals(awg);
                        construction.AWG_e = round(Wires.AWG(awg).AWG_e(validCon));
                        construction.A_b = Wires.AWG(awg).A_b(validCon);
                        construction.D_o = Wires.AWG(awg).D_o(validCon);
                        construction.N_bpl = floor(b./Wires.AWG(awg).D_o(validCon));
                        
                        if isempty(fieldnames(res.constructions))
                            res.constructions.N_s = [];
                            res.constructions.AWG = [];
                            res.constructions.AWG_e = [];
                            res.constructions.A_b = [];
                            res.constructions.D_o = [];
                            res.constructions.N_bpl = [];
                        end
                        
                        res.constructions(validCons) = construction;
                    end
                end
            end
            
            clear isValid
        end
    end
end