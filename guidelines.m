%% Guidelines
% Computes litz wire sizing guidelines.

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
    % surface-to-volume ratio of 4/d for a cylinder (S/V = pi*d*L/(pi/4*d^2*L))
    Jtarget = J_MAX*4; % per unit diameter approximation
    
    % Litz wire strand gauges and areas
    awgz = [33, 36:2:48];
    [~, Az] = AWG2m(awgz);
    
    % available wire constructions (New England Wire)
    load('Wires.mat', 'Wires')
    
    % table from Sullivan paper (2014) - k is a parameter, d is strand diameter
    % truncated to available N.E.W. products
	k = [203, 771, 1.8e3, 4.4e3, 10e3, 24e3, 54e3, 115e3]*1e9;
	d = [0.180, 0.127, 0.101, 0.080, 0.063, 0.050, 0.040, 0.032]*1e-3;
    
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
        subPos = 1:2:2*nwp;
        
        for p = 1:nwp
            curSubPos = subPos(p);
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
        subPos = 2:2:2*nws;
        
        for s = 1:nws
            curSubPos = subPos(s);
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
    
    %% Functions
    
    % computeResultsP
    % Computes results for primary.    
    function res = computeResultsP(W, res)        
        validCons = 0;       % valid construction counter
        N = W.N;
        delta = W.delta_p;
        NPW = 1;
        atLeastOne = false;
        
        % Skin depth provides range of bundles
        dmax = W.d_pMax;
        [~, idx] = find(d(dmax > d), 1);    % index of closest safe fit
        
        % Find suboptimal strands if optimal strand can't be found, up to 2delta
        while isempty(idx) && dmax < 2*delta
            dmax = dmax*2;
            [~, idx] = find(d(dmax > d), 1);
        end
        
        ArangeS = Az(idx:end);              % strand copper areas
        AWGrangeS = awgz(idx:end);          % AWG values
        
        % Winding area bounds bundle size above
        AmaxTot = Ab*N*res.Amin/ACutot;	% fraction of bobbin window
        AmaxOut = (pi/4)*AmaxTot/N;     % max outer round bundle area for fit
        res.Amax = AmaxOut*(0.68125)^2; % max bundle copper area
        res.AWGmax = m2AWG(res.Amax);   % max bundle effective AWG
        
        % Current density bounds bundle size below
        res.AWGmin = m2AWG(res.Amin);   % min bundle effective AWG
        
        while ~atLeastOne
            % base case
            nS = round(k(idx:end)*delta^2*b/N); % integer number of strands
            nS = nS*NPW;                        % parallel wire bundles
            ArangeB = ArangeS.*nS;              % bundle copper areas for nS

            % +25% number of strands
            nShi = ceil(nS*1.25);
            ArangeBhi = ArangeS.*nShi;

            % -25% number of strands
            nSlo = floor(nS*0.75);
            ArangeBlo = ArangeS.*nSlo;
            
            atLeastOne = any(ArangeBlo > res.Amin) && ...
                         any(ArangeB > res.Amin) && ...
                         any(ArangeBhi > res.Amin);
            
            if ~atLeastOne
                NPW = NPW + 1;
            end
        end
        
        % generate plots of fictitious constructions and bounds
        if isequal(nwp, 1) || isequal(p, 1)
            hf = figure;
            hf.Position = [hf.Position(1) - hf.Position(3), ...
                           hf.Position(2) - hf.Position(4), ...
                           hf.Position(3)*2, hf.Position(4)*2];
        end
        
        if nwp > 1
            subplot(max(nwp, nws), 2, curSubPos)
        else
            subplot(1, 2, 1)
        end
        
        plot(AWGrangeS, ArangeB*1e6, 'DisplayName', 'Nominal n_{Strands}')
        hold on
        plot(AWGrangeS, ArangeBhi*1e6, 'DisplayName', '+25%')
        plot(AWGrangeS, ArangeBlo*1e6, 'DisplayName', '-25%')
        hl1 = refline(0, res.Amin*1e6);
        hl1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        text(AWGrangeS(1), res.Amin*1e6, ...
             ['A_{min} = ', num2str(res.Amin*1e6, '%.1f'), ' mm^2', ' (', num2str(round(res.AWGmin)), ' AWG)'], ...
             'VerticalAlignment', 'bottom')
        hl2 = refline(0, res.Amax*1e6);
        hl2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        text(AWGrangeS(1), res.Amax*1e6, ...
             ['A_{max} = ', num2str(res.Amax*1e6, '%.1f'), ' mm^2', ' (', num2str(round(res.AWGmax)), ' AWG)'], ...
             'VerticalAlignment', 'top')
        hold off
        grid on
        
        if nwp > 1
            title(['Area Bounds for Primary ', num2str(p)])
        else
            title('Area Bounds for Primary')
        end
        
        ylabel('Bundle Area [mm^2]')
        
        if isequal(nwp, 1) || isequal(p, nwp)
            xlabel('Strand Gauge [AWG]')
        end
        
        if isequal(nwp, 1) || isequal(p, 1)
            legend('show')
        end
        
        for awg = idx:length(AWGrangeS) + idx - 1
            isValid = NPW*Wires.AWG(awg).A_b >= res.Amin & ...
                      NPW*Wires.AWG(awg).A_b <= res.Amax & ...
                      NPW*Wires.AWG(awg).N_s > nSlo(awg) & ...
                      NPW*Wires.AWG(awg).N_s < nShi(awg);
            
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
                        construction.N_pw = NPW;
                        
                        % prevent fresh constructions from errors
                        if isa(res.constructions, 'double')
                            res.constructions = struct();
                        end
                        
                        if isempty(fieldnames(res.constructions))
                            res.constructions.N_s = [];
                            res.constructions.AWG = [];
                            res.constructions.AWG_e = [];
                            res.constructions.A_b = [];
                            res.constructions.D_o = [];
                            res.constructions.N_bpl = [];
                            res.constructions.N_pw = [];
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
        NPW = 1;
        atLeastOne = false;
        
        % Skin depth provides range of bundles
        dmax = W.d_sMax;
        [~, idx] = find(d(dmax > d), 1);    % index of closest safe fit
        
        % Find suboptimal strands if optimal strand can't be found, up to 2delta
        while isempty(idx) && dmax < 2*delta
            dmax = dmax*2;
            [~, idx] = find(d(dmax > d), 1);
        end
        
        ArangeS = Az(idx:end);              % strand copper areas
        AWGrangeS = awgz(idx:end);          % AWG values
        
        % Winding area bounds bundle size above
        AmaxTot = Ab*N*res.Amin/ACutot;	% fraction of bobbin window
        AmaxOut = (pi/4)*AmaxTot/N;     % max outer round bundle area for fit
        res.Amax = AmaxOut*(0.68125)^2; % max bundle copper area
        res.AWGmax = m2AWG(res.Amax);   % max bundle effective AWG
        
        % Current density bounds bundle size below
        res.AWGmin = m2AWG(res.Amin);   % min bundle effective AWG
        
        while ~atLeastOne
            % base case
            nS = round(k(idx:end)*delta^2*b/N); % integer number of strands
            nS = nS*NPW;                        % parallel wire bundles
            ArangeB = ArangeS.*nS;              % bundle copper areas for nS

            % +25% number of strands
            nShi = ceil(nS*1.25);
            ArangeBhi = ArangeS.*nShi;

            % -25% number of strands
            nSlo = floor(nS*0.75);
            ArangeBlo = ArangeS.*nSlo;
            
            atLeastOne = any(ArangeBlo > res.Amin) && ...
                         any(ArangeB > res.Amin) && ...
                         any(ArangeBhi > res.Amin);
                     
            if ~atLeastOne
                NPW = NPW + 1;
            end
        end
        
        if nws > 1
            subplot(max(nwp, nws), 2, curSubPos)
        else
            subplot(1, 2, 2)
        end
        
        plot(AWGrangeS, ArangeB*1e6, 'DisplayName', 'Nominal n_{Strands}')
        hold on
        plot(AWGrangeS, ArangeBhi*1e6, 'DisplayName', '+25%')
        plot(AWGrangeS, ArangeBlo*1e6, 'DisplayName', '-25%')
        hl1 = refline(0, res.Amin*1e6);
        hl1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        text(AWGrangeS(1), res.Amin*1e6, ...
             ['A_{min} = ', num2str(res.Amin*1e6, '%.1f'), ' mm^2', ' (', num2str(round(res.AWGmin)), ' AWG)'], ...
             'VerticalAlignment', 'bottom')
        hl2 = refline(0, res.Amax*1e6);
        hl2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        text(AWGrangeS(1), res.Amax*1e6, ...
             ['A_{max} = ', num2str(res.Amax*1e6, '%.1f'), ' mm^2', ' (', num2str(round(res.AWGmax)), ' AWG)'], ...
             'VerticalAlignment', 'top')
        hold off
        grid on
        
        if nws > 1
            title(['Area Bounds for Secondary ', num2str(s)])
        else
            title('Area Bounds for Secondary')
        end
        
        ylabel('Bundle Area [mm^2]')
        
        if isequal(nws, 1) || isequal(s, nws)
            xlabel('Strand Gauge [AWG]')
        end
        
        for awg = idx:length(AWGrangeS) + idx - 1
            isValid = NPW*Wires.AWG(awg).A_b >= res.Amin & ...
                      NPW*Wires.AWG(awg).A_b <= res.Amax & ...
                      NPW*Wires.AWG(awg).N_s > nSlo(awg) & ...
                      NPW*Wires.AWG(awg).N_s < nShi(awg);
            
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
                        construction.N_pw = NPW;
                        
                        % prevent fresh constructions from errors
                        if isa(res.constructions, 'double')
                            res.constructions = struct();
                        end
                        
                        if isempty(fieldnames(res.constructions))
                            res.constructions.N_s = [];
                            res.constructions.AWG = [];
                            res.constructions.AWG_e = [];
                            res.constructions.A_b = [];
                            res.constructions.D_o = [];
                            res.constructions.N_bpl = [];
                            res.constructions.N_pw = [];
                        end
                        
                        res.constructions(validCons) = construction;
                    end
                end
            end
            
            clear isValid
        end
    end
end