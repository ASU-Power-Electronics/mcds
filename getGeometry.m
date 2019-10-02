%% getGeometry
% Gets core geometry values using standardized entry format (Ferroxcube).  Works
% for most any core (EPCOS/TDK, Mag. Inc., etc.). Returns structure to add to
% core structure, approximate surface area, and other calculated measures.

function [geometry, As, measures] = getGeometry(name)
    geometry = struct('A', 0, ...
                      'B', 0, ...
                      'C', 0, ...
                      'D', 0, ...
                      'D1', 0, ...
                      'D2', 0, ...
                      'D3', 0, ...
                      'D4', 0, ...
                      'E', 0, ...
                      'F', 0, ...
                      'G', 0, ...
                      'H', 0, ...
                      'H1', 0, ...
                      'H2', 0, ...
                      'H3', 0, ...
                      'ID', 0, ...
                      'OD', 0);
	measures = struct();
                  
    if strcmp(name(1), 'T')
        inputStr = {'H', 'ID', 'OD'};
        defInput = {'0', '0', '0'};
    elseif strcmp(name(1:2), 'RM')
        inputStr = {'A', 'B', 'C', 'D2', 'D3', 'E', 'G', 'H1', 'H2', 'H3'};
        defInput = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};
    elseif strcmp(name(1:2), 'PQ')
        inputStr = {'A', 'B', 'C', 'D2', 'D3', 'E', 'F', 'H1', 'H2'};
        defInput = {'0', '0', '0', '0', '0', '0', '0', '0', '0'};
    elseif strcmp(name(1), 'P') % P or PT
        inputStr = {'A', 'B', 'D1', 'D2', 'D3', 'D4', 'H1', 'H2'};
        defInput = {'0', '0', '0', '0', '0', '0', '0', '0'};
    elseif strcmp(name(1:2), 'ER')
        inputStr = {'A', 'B', 'C', 'D3', 'E', 'F'};
        defInput = {'0', '0', '0', '0', '0', '0'};
    elseif strcmp(name(1:3), 'ETD')
        inputStr = {'A', 'B', 'D2', 'D3', 'E', 'F'};
        defInput = {'0', '0', '0', '0', '0', '0'};
    elseif strcmp(name(1:3), 'EFD')
        inputStr = {'A', 'C', 'E', 'F', 'G', 'H1', 'H2'};
        defInput = {'0', '0', '0', '0', '0', '0', '0'};
    else % E, EC, EPLT
        inputStr = {'A', 'B', 'C', 'D', 'E', 'F'};
        defInput = {'0', '0', '0', '0', '0', '0'};
    end
    
    resp = inputdlg(inputStr, 'Core Measurements', 1, defInput, 'on');
    
    for r = 1:numel(resp)
        geometry.(inputStr{r}) = str2double(resp{r})*1e-3;
    end
    
    if strcmp(name(1), 'T')
        H = geometry.H;
        ID = geometry.ID;
        OD = geometry.OD;
        
        measures.window.height = ID;
        measures.window.breadth = pi*ID;
        measures.d_center = (OD - ID)/2;
        measures.d_center2 = H;
        
        As = pi*((ID + OD)*H + (OD^2 - ID^2)/2);
    elseif strcmp(name(1:2), 'RM')
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D2 = geometry.D2;
        D3 = geometry.D3;
        E = geometry.E;
        G = geometry.G;
        H1 = geometry.H1;
        H2 = geometry.H2;
        H3 = geometry.H3;
        
        measures.window.height = (D2 - D3)/2;
        measures.window.breadth = H2;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = A^2 - B*(sqrt(2)*A - G)/2 - 71*E*(sqrt(2)*A - C)/100 + 2*B*H3 + ...
             pi*D3*H2 + 4*((H1 + H3)/2*(A - (B + 11*E/10)/sqrt(2))) + ...
             2*D2*H2*atan2(sqrt(D2^2 - E^2), E) + pi*D2^2/4 - ...
             9*E*(sqrt(2)*A - C)/25 - D2^2*atan2(E, sqrt(D2^2 - E^2)) + ...
             E*sqrt(D2^2 - E^2)/2;
    elseif strcmp(name(1:2), 'PQ')
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D2 = geometry.D2;
        D3 = geometry.D3;
        E = geometry.E;
        H1 = geometry.H1;
        H2 = geometry.H2;
        
        measures.window.height = (D2 - D3)/2;
        measures.window.breadth = H2;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = 2*A*B - 2*(C + E)*(B - D3) + 2*B*H2 + pi*D3*H1 + 2*A*H2 - ...
             2*C*H1 + D2*H1*atan2(B, C) + (D2^2 - D3^2)*atan2(B, C);
    elseif strcmp(name(1:2), 'PT')
        A = geometry.A;
        B = geometry.B;
        D1 = geometry.D1;
        D3 = geometry.D3;
        D4 = geometry.D4;
        H1 = geometry.H1;
        H2 = geometry.H2;
        
        measures.window.height = (D1 - D3)/2;
        measures.window.breadth = H2;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = pi*(A^2/2 + A*H1)/2 + pi*D4*H1/2 - pi*D4^2/4 + ...
             pi*(D1^2 - D3^2)/4 + pi*(D1 + D3)*H2/2 + B*sqrt(A^2 - B^2)/2 + ...
             A^2*asin(B/A) + B*sqrt(D1^2 - B^2)/2 + D1^2*asin(B/D1) + ...
             A*H1*atan2(B, sqrt(A^2 - B^2)) + ...
             D1*H2*atan2(B, sqrt(D1^2 - B^2)) + ...
             pi*(D4*H1/2 + D3*H2/2 - D4^2/4 - D3^2/4);
    elseif strcmp(name(1), 'P')
        D1 = geometry.D1;
        D2 = geometry.D2;
        D3 = geometry.D3;
        D4 = geometry.D4;
        H1 = geometry.H1;
        H2 = geometry.H2;
        
        measures.window.height = (D2 - D3)/2;
        measures.window.breadth = H2;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = pi*(D1^2/2 + D1*H1) + pi*D4*H1 - pi*D4^2/2 + ...
             pi*(D2^2 - D3^2)/2 + pi*(D2 + D3)*H2;
    elseif strcmp(name(1:2), 'ER')
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D3 = geometry.D3;
        E = geometry.E;
        F = geometry.F;
        
        measures.window.height = (C - D3)/2;
        measures.window.breadth = 2*F;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = 2*A*B + 4*B*(E + F) + 2*pi*D3*F + 4*(A*E - F*C) + 2*B*C - pi*D3^2/2;
    elseif strcmp(name(1:3), 'ETD')
        A = geometry.A;
        B = geometry.B;
        D2 = geometry.D2;
        D3 = geometry.D3;
        E = geometry.E;
        F = geometry.F;
        
        measures.window.height = (D2 - D3)/2;
        measures.window.breadth = 2*F;
        measures.d_center = D3;
        measures.d_center2 = D3;
        
        As = 2*A*B + 4*B*E + 4*A*(E - F) + ...
             8*F*(A/2 - D2*cos(2*atan2(B, D2))/2) + 2*pi*D3*F + ...
             4*F*D2*atan2(B, D2) + B*sqrt(D2^2 - B^2) + 2*D2^2*asin(B/D2) - ...
             pi*D3^2/2;
    elseif strcmp(name(1:4), 'EPLT')
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D = geometry.D;
        E = geometry.E;
        F = geometry.F;
        
        measures.window.height = (B - C)/2;
        measures.window.breadth = 2*E;
        measures.d_center = C;
        measures.d_center2 = F;
        
        As = 2*A*(D + E + F) + 2*E*(C - B) + 2*D*F + 6*E*F + 2*F*(B - C);
    elseif strcmp(name(1:3), 'EFD')
        A = geometry.A;
        C = geometry.C;
        E = geometry.E;
        F = geometry.F;
        G = geometry.G;
        H1 = geometry.H1;
        H2 = geometry.H2;
        
        measures.window.height = (E - F)/2;
        measures.window.breadth = 2*H2;
        measures.d_center = F;
        measures.d_center2 = G;
        
        As = 4*C*(H1 + H2) + 4*H2*(A - E) + 4*(A + C - G)*(H1 - H2) + 4*H2*(F + G) + 2*C*(E - F);
    elseif strcmp(name(1:2), 'EC')
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D = geometry.D;
        E = geometry.E;
        F = geometry.F;
        
        measures.window.height = (B - C)/2;
        measures.window.breadth = 2*E;
        measures.d_center = C;
        measures.d_center2 = C;
        
        As = 2*A*F + 4*F*D + 4*A*(D - E) + ...
             8*E*(A/2 - B*cos(2*atan2(F, B))/2) + 2*pi*C*E + ...
             4*E*B*atan2(F, B) + F*sqrt(B^2 - F^2) + 2*B^2*asin(F/B) - ...
             pi*C^2/2;
    else % E
        A = geometry.A;
        B = geometry.B;
        C = geometry.C;
        D = geometry.D;
        E = geometry.E;
        F = geometry.F;
        
        measures.window.height = (B - C)/2;
        measures.window.breadth = 2*E;
        measures.d_center = C;
        measures.d_center2 = F;
        
        As = 2*F*(A + 2*D + 4*E) + 4*E*(A - B + C) + 4*A*(D - E);
    end
end