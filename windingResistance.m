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

    FlitzL = 1.05; % twisting pitch correction for litz wire
    eta = 0.75;  % porosity factor (approximation ~5% for small d_s)
    
    % Strand outer diameters from Sullivan (1999)
    alf = 1.12;
    bet = 0.97;
    [dref, ~] = AWG2m(40);
    
    P = 0; % initialize total loss in all windings to 0 for summation
    
    % create inputdlg structure for iterated input dialog
    inputStruct = struct();
    inputStruct.prompt = {'Primary winding', 'termination (flag) length [mm]:'; ...
                          'Secondary winding', 'termination (flag) length [mm]:'};
    inputStruct.dlgtitle = 'Flag Lengths';
    inputStruct.definput = 100;
    inputStruct.np = nwp;
    inputStruct.ns = nws;
    
    % gather responses for flag lengths
    answer = createWindingInput(inputStruct);
    
    % Primary winding(s)
    if isequal(nwp, 1)
        thisP.flagLength = str2double(answer{1})*1e-3;
        
        % geometry and counts
        ds = thisP.d_s;
        dl = dref*alf*(thisP.d_s/dref)^bet; % layer diameter (strand outer diameter)
        k = thisP.N_s/thisP.NPW; % number of strands per bundle
        length = thisP.length + thisP.flagLength; % length of winding
        Fi = thisP.F_i; % interleaving factor, correction to Dowell

        % DC resistance and loss
        RsDC = 4*RHO_CU*length*FlitzL/(pi*ds^2); % DC resistance of strand
        RwDC = RsDC/(k*thisP.NPW); % DC resistance of winding
        Pw = RwDC*thisP.I_pRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisP.N_L; % number of layers
        deltaw = thisP.delta_p; % skin depth of current in winding
        Nll = Nl*sqrt(k)/Fi; % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR_d = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)));
        FR_p = A*((2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));
        FR = FR_d + FR_p;
        
        thisP.FR_d = FR_d;
        thisP.FR_p = FR_p;
        thisP.FR = FR;

        % resistance and copper loss
        thisP.R_DC = RwDC;
        thisP.R = RwDC*FR;
        thisP.P_Cu = Pw*FR;
        P = P + thisP.P_Cu; % sum to total copper loss
    else
        for p = 1:nwp
            thisP(p).flagLength = str2double(answer{p})*1e-3;
            
            % geometry and counts
            ds = thisP(p).d_s;
            dl = dref*alf*(thisP(p).d_s/dref)^bet; % layer diameter (strand outer diameter)
            k = thisP(p).N_s/thisP(p).NPW; % number of strands per bundle
            length = thisP(p).length + thisP(p).flagLength; % length of winding
            Fi = thisP(p).F_i; % interleaving factor, correction to Dowell

            % DC resistance and loss
            RsDC = 4*RHO_CU*length*FlitzL/(pi*ds^2); % DC resistance of strand
            RwDC = RsDC/(k*thisP(p).NPW); % DC resistance of winding
            Pw = RwDC*thisP(p).I_pRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisP(p).N_L; % number of layers
            deltaw = thisP(p).delta_p; % skin depth of current in winding
            Nll = Nl*sqrt(k)/Fi; % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR_d = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)));
            FR_p = A*((2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));
            FR = FR_d + FR_p;

            thisP(p).FR_d = FR_d;
            thisP(p).FR_p = FR_p;
            thisP(p).FR = FR;

            thisP(p).R_DC = RwDC;
            thisP(p).R = RwDC*FR;
            thisP(p).P_Cu = Pw*FR;
            P = P + thisP(p).P_Cu;
        end
    end
    
    % Secondary winding(s)
    if isequal(nws, 1)
        thisS.flagLength = str2double(answer{2})*1e-3;
        
        % geometry and counts
        ds = thisS.d_s;
        dl = dref*alf*(thisS.d_s/dref)^bet; % layer diameter (strand outer diameter)
        k = thisS.N_s/thisS.NPW; % number of strands per bundle
        length = thisS.length + thisS.flagLength; % length of winding
        Fi = thisS.F_i; % interleaving factor, correction to Dowell

        % DC resistance and loss
        RsDC = 4*RHO_CU*length*FlitzL/(pi*ds^2); % DC resistance of strand
        RwDC = RsDC/(k*thisS.NPW); % DC resistance of winding
        Pw = RwDC*thisS.I_sRMS^2; % DC loss of winding

        % Dowell's equation
        Nl = thisS.N_L; % number of layers
        deltaw = thisS.delta_s; % skin depth of current in winding
        Nll = Nl*sqrt(k)/Fi; % effective number of strand layers (square)

        A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
        FR_d = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)));
        FR_p = A*((2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));
        FR = FR_d + FR_p;
        
        thisS.FR_d = FR_d;
        thisS.FR_p = FR_p;
        thisS.FR = FR;

        % resistance and copper loss
        thisS.R_DC = RwDC;
        thisS.R = RwDC*FR;
        thisS.P_Cu = Pw*FR;
        P = P + thisS.P_Cu; % sum to total copper loss
    else
        for s = 1:nws
            thisS(s).flagLength = str2double(answer{p + s})*1e-3;
            
            % geometry and counts
            ds = thisS(s).d_s;
            dl = dref*alf*(thisS(s).d_s/dref)^bet; % layer diameter (strand outer diameter)
            k = thisS(s).N_s/thisS(s).NPW; % number of strands per bundle
            length = thisS(s).length + thisS(s).flagLength; % length of winding
            Fi = thisS(s).F_i; % interleaving factor, correction to Dowell

            % DC resistance and loss
            RsDC = 4*RHO_CU*length*FlitzL/(pi*ds^2); % DC resistance of strand
            RwDC = RsDC/(k*thisS(s).NPW); % DC resistance of winding
            Pw = RwDC*thisS(s).I_sRMS^2; % DC loss of winding

            % Dowell's equation
            Nl = thisS(s).N_L; % number of layers
            deltaw = thisS(s).delta_s; % skin depth of current in winding
            Nll = Nl*sqrt(k)/Fi; % effective number of strand layers (square)

            A = (pi/4)^0.75*dl/deltaw*sqrt(eta);
            FR_d = A*((sinh(2*A) + sin(2*A))/(cosh(2*A) - cos(2*A)));
            FR_p = A*((2*(Nll^2 - 1)/3)*((sinh(A) - sin(A))/(cosh(A) + cos(A))));
            FR = FR_d + FR_p;

            thisS(s).FR_d = FR_d;
            thisS(s).FR_p = FR_p;
            thisS(s).FR = FR;

            thisS(s).R_DC = RwDC;
            thisS(s).R = RwDC*FR;
            thisS(s).P_Cu = Pw*FR;
            P = P + thisS(s).P_Cu;
        end
    end
    
    Winding.primary = orderfields(thisP);
    Winding.secondary = orderfields(thisS);
end