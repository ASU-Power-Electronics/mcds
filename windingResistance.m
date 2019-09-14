%% windingResistance
% Approximation of winding resistance and loss for litz wire-wound transformers.
% Based on formulation in Kazimierczuk.

%TODO: add functionality for varying wire types (foil, planar, etc.)
%TODO: calculate eta instead of declaring it

function [Winding, P] = windingResistance(Winding)
    % variable assignments for brevity
    global RHO_CU
    
    thisP = Winding.primary;
    thisS = Winding.secondary;
    [~, nwp] = size(thisP);
    [~, nws] = size(thisS);

    eta = 0.75;  % porosity factor (approximation ~5% for small d_s)
    
    % Strand outer diameters from Sullivan (1999)
    alf = 1.12;
    bet = 0.97;
    [dref, ~] = AWG2m(40);
    
    P = 0; % initialize total loss in all windings to 0 for summation
    
    % Primary winding(s)
    if isequal(nwp, 1)
        % request winding termination length
        answer = inputdlg('Primary winding termination (flag) length [mm]:', ...
                          'Primary Flag Length', ...
                          1, ...
                          {'100'});
        thisP.flagLength = str2double(answer{1})*1e-3;
        
        % geometry and counts
        ds = thisP.d_s;
        dl = dref*alf*(thisP.d_s/dref)^bet; % layer diameter (strand outer diameter)
        k = thisP.N_s; % number of strands
        length = thisP.length + thisP.flagLength; % length of winding

        % DC resistance and loss
        RsDC = 4*RHO_CU*length/(pi*ds^2); % DC resistance of strand
        RwDC = RsDC/(k*thisP.NPW); % DC resistance of winding
        Pw = RwDC*thisP.I_pRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisP.N_L; % number of layers
        deltaw = thisP.delta_p; % skin depth of current in winding
        Nll = Nl*sqrt(k); % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

        % resistance and copper loss
        thisP.R_DC = RwDC;
        thisP.R = RwDC*FR;
        thisP.P_Cu = Pw*FR;
        P = P + thisP.P_Cu; % sum to total copper loss
    else
        for w = 1:nwp
            % request winding termination length
            answer = inputdlg(sprintf('Primary winding %d termination (flag) length [mm]:', w), ...
                              sprintf('Primary %d Flag Length', w), ...
                              1, ...
                              {'100'});
            thisP(w).flagLength = str2double(answer{1})*1e-3;
            
            % geometry and counts
            ds = thisP(w).d_s;
            dl = dref*alf*(thisP(w).d_s/dref)^bet; % layer diameter (strand outer diameter)
            k = thisP(w).N_s; % number of strands
            length = thisP(w).length + thisP(w).flagLength; % length of winding

            % DC resistance and loss
            RsDC = 4*RHO_CU*length/(pi*ds^2); % DC resistance of strand
            RwDC = RsDC/(k*thisP(w).NPW); % DC resistance of winding
            Pw = RwDC*thisP(w).I_pRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisP(w).N_L; % number of layers
            deltaw = thisP(w).delta_p; % skin depth of current in winding
            Nll = Nl*sqrt(k); % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

            thisP(w).R_DC = RwDC;
            thisP(w).R = RwDC*FR;
            thisP(w).P_Cu = Pw*FR;
            P = P + thisP(w).P_Cu;
        end
    end
    
    % Secondary winding(s)
    if isequal(nws, 1)
        % request winding termination length
        answer = inputdlg('Secondary winding termination (flag) length [mm]:', ...
                          'Secondary Flag Length', ...
                          1, ...
                          {'100'});
        thisS.flagLength = str2double(answer{1})*1e-3;
        
        % geometry and counts
        ds = thisS.d_s;
        dl = dref*alf*(thisS.d_s/dref)^bet; % layer diameter (strand outer diameter)
        k = thisS.N_s; % number of strands
        length = thisS.length + thisS.flagLength; % length of winding

        % DC resistance and loss
        RsDC = 4*RHO_CU*length/(pi*ds^2); % DC resistance of strand
        RwDC = RsDC/(k*thisS.NPW); % DC resistance of winding
        Pw = RwDC*thisS.I_sRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisS.N_L; % number of layers
        deltaw = thisS.delta_s; % skin depth of current in winding
        Nll = Nl*sqrt(k); % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

        % resistance and copper loss
        thisS.R_DC = RwDC;
        thisS.R = RwDC*FR;
        thisS.P_Cu = Pw*FR;
        P = P + thisS.P_Cu; % sum to total copper loss
    else
        for w = 1:nws
            % request winding termination length
            answer = inputdlg(sprintf('Secondary winding %d termination (flag) length [mm]:', w), ...
                              sprintf('Secondary %d Flag Length', w), ...
                              1, ...
                              {'100'});
            thisS(w).flagLength = str2double(answer{1})*1e-3;
            
            % geometry and counts
            ds = thisS(w).d_s;
            dl = dref*alf*(thisS(w).d_s/dref)^bet; % layer diameter (strand outer diameter)
            k = thisS(w).N_s; % number of strands
            length = thisS(w).length + thisS(w).flagLength; % length of winding

            % DC resistance and loss
            RsDC = 4*RHO_CU*length/(pi*ds^2); % DC resistance of strand
            RwDC = RsDC/(k*thisS(w).NPW); % DC resistance of winding
            Pw = RwDC*thisS(w).I_sRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisS(w).N_L; % number of layers
            deltaw = thisS(w).delta_s; % skin depth of current in winding
            Nll = Nl*sqrt(k); % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

            thisS(w).R_DC = RwDC;
            thisS(w).R = RwDC*FR;
            thisS(w).P_Cu = Pw*FR;
            P = P + thisS(w).P_Cu;
        end
    end
    
    Winding.primary = orderfields(thisP);
    Winding.secondary = orderfields(thisS);
end