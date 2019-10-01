%% arrangeWindings
% Arranges windings in window.  Allows interleaving of windings, but only full
% interleaving (PSPS...PS/SP), since it is most advantageous.

function [thisP, thisS] = arrangeWindings(thisP, nwp, thisS, nws, Wgap, Core)
    N_Lp = 0;
    N_Ls = 0;
    
    if nwp > 1
        for idx = 1:nwp
            if isequal(thisP(idx).bifilar, 2)
                N_Lp = N_Lp + thisP(idx).N_L/2;
            else
                N_Lp = N_Lp + thisP(idx).N_L;
            end
        end
    else
        N_Lp = N_Lp + thisP.N_L;
    end
    
    if nws > 1
        for idx = 1:nws
            if isequal(thisS(idx).bifilar, 2)
                N_Ls = N_Ls + thisS(idx).N_L/2;
            else
                N_Ls = N_Ls + thisS(idx).N_L;
            end
        end
    else
        N_Ls = N_Ls + thisS.N_L;
    end

    % check conditions for interleaved winding eligibility
    isEligibleForInterleave = (isequal(N_Lp, N_Ls) || ...
                               isequal(N_Lp, N_Ls + 1)) && ...
                               N_Lp > 1;

    if isEligibleForInterleave
        windingStructure = questdlg('Simple or interleaved winding structure?', ...
                                    'Winding Structure Selection', ...
                                    'Simple', 'Interleaved', 'Simple');
    else
        windingStructure = 'Simple';
    end

    % check for round/rectangular cross-section and assign approx. initial diameter
    if strcmp(Core.c_leg_type, 'round')
        lastDiameter = Core.d_center + 2*(Core.window.height - Core.bobbin.height);
    else % rectangular
        lastDiameter = sqrt(Core.d_center^2 + Core.d_center2^2) + 2*(Core.window.height - Core.bobbin.height);
        theta = atan2(Core.d_center2, Core.d_center); % for similar triangle later
    end

    switch windingStructure
        case 'Simple' % direct concentric (PPP...SSS...)
            [thisP(:).F_i] = 1; % interleaving factor, correction to Dowell
            [thisS(:).F_i] = 1;
            
            if nwp > 1
                bifMatch = 0;

                for idx = 1:nwp
                    if isequal(thisP(idx).bifilar, 2)
                        if bifMatch
                            bifMatch = 0;
                            thisP(idx).diameter = lastDiameter - thisP(idx - 1).d_o*thisP(idx - 1).N_L;
                            
                            % length of winding
                            if strcmp(Core.c_leg_type, 'rectangular') % Pythagoras
                                c = thisP(idx).diameter;
                                a = c*cos(theta);
                                b = c*sin(theta);
                                thisP(idx).length = 2*(a + b)*thisP(idx).N;
                            else
                                thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                            end
                        else
                            bifMatch = 1;
                            % diameter at center of winding (for use in L calculations)
                            thisP(idx).diameter = lastDiameter + thisP(idx).d_o*thisP(idx).N_L;
                            % outer diameter of winding (for next winding)
                            lastDiameter = thisP(idx).diameter + thisP(idx).d_o*thisP(idx).N_L;
                            % inter-winding space
                            lastDiameter = lastDiameter + Wgap;
                            
                            if strcmp(Core.c_leg_type, 'rectangular')
                                c = thisP(idx).diameter;
                                a = c*cos(theta);
                                b = c*sin(theta);
                                thisP(idx).length = 2*(a + b)*thisP(idx).N;
                            else
                                thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                            end
                        end
                    else
                        thisP(idx).diameter = lastDiameter + thisP(idx).d_o*thisP(idx).N_L;
                        lastDiameter = thisP(idx).diameter + thisP(idx).d_o*thisP(idx).N_L;
                        lastDiameter = lastDiameter + Wgap;
                        
                        if strcmp(Core.c_leg_type, 'rectangular')
                            c = thisP(idx).diameter;
                            a = c*cos(theta);
                            b = c*sin(theta);
                            thisP(idx).length = 2*(a + b)*thisP(idx).N;
                        else
                            thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                        end
                    end
                end
            else % single primary
                thisP.diameter = lastDiameter + thisP.d_o*thisP.N_L;
                lastDiameter = thisP.diameter + thisP.d_o*thisP.N_L;
                lastDiameter = lastDiameter + Wgap;
                
                if strcmp(Core.c_leg_type, 'rectangular')
                    c = thisP.diameter;
                    a = c*cos(theta);
                    b = c*sin(theta);
                    thisP.length = 2*(a + b)*thisP.N;
                else
                    thisP.length = pi*thisP.diameter*thisP.N;
                end
            end

            if nws > 1
                bifMatch = 0;

                for idx = 1:nws
                    if isequal(thisS(idx).bifilar, 2)
                        if bifMatch
                            bifMatch = 0;
                            thisS(idx).diameter = lastDiameter - thisS(idx - 1).d_o*thisS(idx - 1).N_L;
                
                            if strcmp(Core.c_leg_type, 'rectangular')
                                c = thisS(idx).diameter;
                                a = c*cos(theta);
                                b = c*sin(theta);
                                thisS(idx).length = 2*(a + b)*thisS(idx).N;
                            else
                                thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                            end
                        else
                            bifMatch = 1;
                            thisS(idx).diameter = lastDiameter + thisS(idx).d_o*thisS(idx).N_L;
                            lastDiameter = thisS(idx).diameter + thisS(idx).d_o*thisS(idx).N_L;
                            
                            if idx < nws
                                lastDiameter = lastDiameter + Wgap;
                            end
                
                            if strcmp(Core.c_leg_type, 'rectangular')
                                c = thisS(idx).diameter;
                                a = c*cos(theta);
                                b = c*sin(theta);
                                thisS(idx).length = 2*(a + b)*thisS(idx).N;
                            else
                                thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                            end
                        end
                    else
                        thisS(idx).diameter = lastDiameter + thisS(idx).d_o*thisS(idx).N_L;
                        lastDiameter = thisS(idx).diameter + thisS(idx).d_o*thisS(idx).N_L;
                        
                        if idx < nws
                            lastDiameter = lastDiameter + Wgap;
                        end
                
                        if strcmp(Core.c_leg_type, 'rectangular')
                            c = thisS(idx).diameter;
                            a = c*cos(theta);
                            b = c*sin(theta);
                            thisS(idx).length = 2*(a + b)*thisS(idx).N;
                        else
                            thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                        end
                    end
                end
            else
                thisS.diameter = lastDiameter + thisS.d_o*thisS.N_L;
                lastDiameter = thisS.diameter + thisS.d_o*thisS.N_L;
                
                if strcmp(Core.c_leg_type, 'rectangular')
                    c = thisS.diameter;
                    a = c*cos(theta);
                    b = c*sin(theta);
                    thisS.length = 2*(a + b)*thisS.N;
                else
                    thisS.length = pi*thisS.diameter*thisS.N;
                end
            end
            
            % quick check that winding height < window height
            if strcmp(Core.c_leg_type, 'round')
                if ~(lastDiameter < Core.d_center + 2*Core.window.height)
                    warning('Windings exceed window height as arranged.')
                end
            else
                hWp = [thisP(:).N_L]'*[thisP(:).d_o];
                hWs = [thisS(:).N_L]'*[thisS(:).d_o];
                
                if ~(hWp + hWs < Core.window.height)
                    warning('Windings exceed window height as arranged.')
                end
            end
        case 'Interleaved' % alternating concentric (PSPS...PS/SP)
            lastPS = 'S'; % always start with primary (higher power due to loss)
            pIDX = 1;
            pCurrentLayer = 1;
            sIDX = 1;
            sCurrentLayer = 1;
            
             % interleaving factor, correction to Dowell
            for p = 1:nwp
                thisP(p).F_i = N_Lp;
            end
            
            for s = 1:nws
                thisS(s).F_i = N_Ls;
            end
            
            % build layer diameter array
            for layer = 1:(N_Lp + N_Ls)
                switch lastPS
                    case 'S' % primary layer
                        if nwp > 1
                            if isequal(pCurrentLayer, 1)
                                d_Pin = lastDiameter;
                            end
                            
                            if pCurrentLayer < thisP(pIDX).N_L                                
                                lastDiameter = lastDiameter + 2*thisP(pIDX).d_o;
                                pCurrentLayer = pCurrentLayer + 1;
                            else
                                d_Pout = lastDiameter + 2*thisP(pIDX).d_o;
                                lastDiameter = d_Pout + Wgap;
                                thisP(pIDX).diameter = (d_Pin + d_Pout)/2;
                
                                if strcmp(Core.c_leg_type, 'rectangular')
                                    c = thisP(pIDX).diameter;
                                    a = c*cos(theta);
                                    b = c*sin(theta);
                                    thisP(pIDX).length = 2*(a + b)*thisP(pIDX).N;
                                else
                                    thisP(pIDX).length = pi*thisP(pIDX).diameter*thisP(pIDX).N;
                                end
                                
                                if isequal(thisP(pIDX).bifilar, 2)
                                    thisP(pIDX + 1).diameter = thisP(pIDX).diameter;
                                    thisP(pIDX + 1).length = thisP(pIDX).length;
                                    pIDX = pIDX + 2;
                                else
                                    pIDX = pIDX + 1;
                                end
                                
                                pCurrentLayer = 1;
                            end
                        else
                            if isequal(pCurrentLayer, 1)
                                d_Pin = lastDiameter;
                            end
                            
                            if pCurrentLayer < thisP.N_L                                
                                lastDiameter = lastDiameter + 2*thisP.d_o;
                                pCurrentLayer = pCurrentLayer + 1;
                            else
                                d_Pout = lastDiameter + 2*thisP.d_o;
                                lastDiameter = d_Pout + Wgap;
                                thisP.diameter = (d_Pin + d_Pout)/2;
                
                                if strcmp(Core.c_leg_type, 'rectangular')
                                    c = thisP.diameter;
                                    a = c*cos(theta);
                                    b = c*sin(theta);
                                    thisP.length = 2*(a + b)*thisP.N;
                                else
                                    thisP.length = pi*thisP.diameter*thisP.N;
                                end
                            end
                        end
                        
                        lastPS = 'P';
                    case 'P' % secondary layer
                        if nws > 1
                            if isequal(sCurrentLayer, 1)
                                d_Sin = lastDiameter;
                            end
                            
                            if sCurrentLayer < thisS(sIDX).N_L                                
                                lastDiameter = lastDiameter + 2*thisS(sIDX).d_o;
                                sCurrentLayer = sCurrentLayer + 1;
                            else
                                d_Sout = lastDiameter + 2*thisS(sIDX).d_o;
                                lastDiameter = d_Sout + Wgap;
                                thisS(sIDX).diameter = (d_Sin + d_Sout)/2;
                                thisS(sIDX).length = pi*thisS(sIDX).diameter*thisS(sIDX).N;
                
                                if strcmp(Core.c_leg_type, 'rectangular')
                                    c = thisS(sIDX).diameter;
                                    a = c*cos(theta);
                                    b = c*sin(theta);
                                    thisS(sIDX).length = 2*(a + b)*thisS(sIDX).N;
                                else
                                    thisS(sIDX).length = pi*thisS(sIDX).diameter*thisS(sIDX).N;
                                end
                                
                                if isequal(thisS(sIDX).bifilar, 2)
                                    thisS(sIDX + 1).diameter = thisS(sIDX).diameter;
                                    thisS(sIDX + 1).length = thisS(sIDX).length;
                                    sIDX = sIDX + 2;
                                else
                                    sIDX = sIDX + 1;
                                end
                                
                                sCurrentLayer = 1;
                            end
                        else
                            if isequal(sCurrentLayer, 1)
                                d_Sin = lastDiameter;
                            end
                            
                            if sCurrentLayer < thisS.N_L                                
                                lastDiameter = lastDiameter + 2*thisS.d_o;
                                sCurrentLayer = sCurrentLayer + 1;
                            else
                                d_Sout = lastDiameter + 2*thisS.d_o;
                                lastDiameter = d_Sout + Wgap;
                                thisS.diameter = (d_Sin + d_Sout)/2;
                
                                if strcmp(Core.c_leg_type, 'rectangular')
                                    c = thisS.diameter;
                                    a = c*cos(theta);
                                    b = c*sin(theta);
                                    thisS.length = 2*(a + b)*thisS.N;
                                else
                                    thisS.length = pi*thisS.diameter*thisS.N;
                                end
                            end
                        end
                        
                        lastPS = 'S';
                end
            end
            
            % quick check that winding height < window height
            if strcmp(Core.c_leg_type, 'round')
                if lastDiameter > Core.d_center + 2*Core.window.height
                    warning('Windings exceed window height as arranged.')
                end
            else
                hWp = [thisP(:).N_L]'*[thisP(:).d_o];
                hWs = [thisS(:).N_L]'*[thisS(:).d_o];
                
                if hWp + hWs > Core.window.height
                    warning('Windings exceed window height as arranged.')
                end
            end
    end
end