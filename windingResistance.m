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
    
    P = 0; % initialize total loss in all windings to 0 for summation
    
    % Primary winding(s)
    if nwp == 1
        % request winding termination length
        answer = inputdlg('Primary winding termination (flag) length [mm]:', ...
                          'Primary Flag Length', ...
                          1, ...
                          {'100'});
        thisP.flagLength = str2double(answer{1})*1e-3;
        
        % geometry and counts
        dl = thisP.d_s; % layer diameter (strand diameter)
        k = thisP.N_s; % number of strands
        length = thisP.length + thisP.flagLength; % length of winding

        % DC resistance and loss
        RsDC = 4*RHO_CU*length/(pi*dl^2); % DC resistance of strand
        RwDC = RsDC/(k*thisP.NPW); % DC resistance of winding
        Pw = RwDC*thisP.I_pRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisP.N_L; % number of layers
        deltaw = thisP.delta_p; % skin depth of current in winding
        Nll = Nl*sqrt(k); % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

        % resistance and copper loss
        thisP.R = RwDC*FR;
        thisP.P_Cu = Pw*FR;
        P = P + thisP.P_Cu; % sum to total copper loss
    else
        for i = 1:nwp
            % request winding termination length
            answer = inputdlg(sprintf('Primary winding %d termination (flag) length [mm]:', i), ...
                              sprintf('Primary %d Flag Length', i), ...
                              1, ...
                              {'100'});
            thisP(i).flagLength = str2double(answer{1})*1e-3;
            
            % geometry and counts
            dl = thisP(i).d_s; % layer diameter (strand diameter)
            k = thisP(i).N_s; % number of strands
            length = thisP(i).length + thisP(i).flagLength; % length of winding

            % DC resistance and loss
            RsDC = 4*RHO_CU*length/(pi*dl^2); % DC resistance of strand
            RwDC = RsDC/(k*thisP(i).NPW); % DC resistance of winding
            Pw = RwDC*thisP(i).I_pRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisP(i).N_L; % number of layers
            deltaw = thisP(i).delta_p; % skin depth of current in winding
            Nll = Nl*sqrt(k); % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

            thisP(i).R = RwDC*FR;
            thisP(i).P_Cu = Pw*FR;
            P = P + thisP(i).P_Cu;
        end
    end
    
    % Secondary winding(s)
    if nws == 1
        % request winding termination length
        answer = inputdlg('Secondary winding termination (flag) length [mm]:', ...
                          'Secondary Flag Length', ...
                          1, ...
                          {'100'});
        thisS.flagLength = str2double(answer{1})*1e-3;
        
        % geometry and counts
        dl = thisS.d_s; % layer diameter (strand diameter)
        k = thisS.N_s; % number of strands
        length = thisS.length + thisS.flagLength; % length of winding

        % DC resistance and loss
        RsDC = 4*RHO_CU*length/(pi*dl^2); % DC resistance of strand
        RwDC = RsDC/(k*thisS.NPW); % DC resistance of winding
        Pw = RwDC*thisS.I_sRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisS.N_L; % number of layers
        deltaw = thisS.delta_s; % skin depth of current in winding
        Nll = Nl*sqrt(k); % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

        % resistance and copper loss
        thisS.R = RwDC*FR;
        thisS.P_Cu = Pw*FR;
        P = P + thisS.P_Cu; % sum to total copper loss
    else
        for i = 1:nws
            % request winding termination length
            answer = inputdlg(sprintf('Secondary winding %d termination (flag) length [mm]:', i), ...
                              sprintf('Secondary %d Flag Length', i), ...
                              1, ...
                              {'100'});
            thisS(i).flagLength = str2double(answer{1})*1e-3;
            
            % geometry and counts
            dl = thisS(i).d_s; % layer diameter (strand diameter)
            k = thisS(i).N_s; % number of strands
            length = thisS(i).length + thisS(i).flagLength; % length of winding

            % DC resistance and loss
            RsDC = 4*RHO_CU*length/(pi*dl^2); % DC resistance of strand
            RwDC = RsDC/(k*thisS(i).NPW); % DC resistance of winding
            Pw = RwDC*thisS(i).I_sRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisS(i).N_L; % number of layers
            deltaw = thisS(i).delta_s; % skin depth of current in winding
            Nll = Nl*sqrt(k); % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)) + (2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));

            thisS(i).R = RwDC*FR;
            thisS(i).P_Cu = Pw*FR;
            P = P + thisS(i).P_Cu;
        end
    end
    
    Winding.primary = orderfields(thisP);
    Winding.secondary = orderfields(thisS);
end