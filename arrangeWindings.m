%% arrangeWindings
% Arranges windings in window.  Allows interleaving of windings, but only full
% interleaving (PSPS...PS/SP), since it is most advantageous.

function [thisP, thisS] = arrangeWindings(thisP, nwp, thisS, nws, Core)
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
    isEligibleForInterleave = isequal(N_Lp, N_Ls) || isequal(N_Lp, N_Ls + 1);

    if isEligibleForInterleave
        windingStructure = questdlg('Simple or interleaved winding structure?', ...
                                    'Winding Structure Selection', ...
                                    'Simple', 'Interleaved', 'Simple');
    else
        windingStructure = 'Simple';
    end

    % check for round/rectangular cross-section and assign approx. initial diameter
    if isequal(Core.d_center, Core.d_center2)
        lastDiameter = Core.d_center + 2*(Core.window.height - Core.bobbin.height);
    else % rectangular
        lastDiameter = 2*sqrt(Core.d_center*Core.d_center2/pi) + 2*(Core.window.height - Core.bobbin.height);
    end

    switch windingStructure
        case 'Simple' % direct concentric (PPP...SSS...)
            if nwp > 1
                bifMatch = 0;

                for idx = 1:nwp
                    if isequal(thisP(idx).bifilar, 2)
                        if bifMatch
                            bifMatch = 0;
                            thisP(idx).diameter = lastDiameter - thisP(idx - 1).d_o*thisP(idx - 1).N_L;
                            thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                        else
                            bifMatch = 1;
                            % diameter at center of winding (for use in L calculations)
                            thisP(idx).diameter = lastDiameter + thisP(idx).d_o*thisP(idx).N_L;
                            % outer diameter of winding (for next winding)
                            lastDiameter = thisP(idx).diameter + thisP(idx).d_o*thisP(idx).N_L;
                            % length using MLT*N
                            thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                        end
                    else
                        thisP(idx).diameter = lastDiameter + thisP(idx).d_o*thisP(idx).N_L;
                        lastDiameter = thisP(idx).diameter + thisP(idx).d_o*thisP(idx).N_L;
                        thisP(idx).length = pi*thisP(idx).diameter*thisP(idx).N;
                    end
                end
            else
                thisP.diameter = lastDiameter + thisP.d_o*thisP.N_L;
                lastDiameter = thisP.diameter + thisP.d_o*thisP.N_L;
                thisP.length = pi*thisP.diameter*thisP.N;
            end

            if nws > 1
                bifMatch = 0;

                for idx = 1:nws
                    if isequal(thisS(idx).bifilar, 2)
                        if bifMatch
                            bifMatch = 0;
                            thisS(idx).diameter = lastDiameter - thisS(idx - 1).d_o*thisS(idx - 1).N_L;
                            thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                        else
                            bifMatch = 1;
                            thisS(idx).diameter = lastDiameter + thisS(idx).d_o*thisS(idx).N_L;
                            lastDiameter = thisS(idx).diameter + thisS(idx).d_o*thisS(idx).N_L;
                            thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                        end
                    else
                        thisS(idx).diameter = lastDiameter + thisS(idx).d_o*thisS(idx).N_L;
                        lastDiameter = thisS(idx).diameter + thisS(idx).d_o*thisS(idx).N_L;
                        thisS(idx).length = pi*thisS(idx).diameter*thisS(idx).N;
                    end
                end
            else
                thisS.diameter = lastDiameter + thisS.d_o*thisS.N_L;
                lastDiameter = thisS.diameter + thisS.d_o*thisS.N_L;
                thisS.length = pi*thisS.diameter*thisS.N;
            end
            
            % quick check that winding height < window height
            if ~(lastDiameter < Core.d_center + 2*Core.window.height)
                warning('Windings exceed window height as arranged.')
            end
        case 'Interleaved' % alternating concentric (PSPS...PS/SP)
            lastPS = 'S'; % always start with primary (higher power due to loss)
            pIDX = 1;
            pCurrentLayer = 1;
            sIDX = 1;
            sCurrentLayer = 1;
            
            % build layer diameter array
            for layer = 1:(N_Lp + N_Ls)
                switch lastPS
                    case 'S' % primary layer
                        if nwp > 1
                            if pCurrentLayer < thisP(pIDX).N_L
                                if isequal(pCurrentLayer, 1)
                                    d_Pin = lastDiameter;
                                end
                                
                                lastDiameter = lastDiameter + 2*thisP(pIDX).d_o;
                            else
                                d_Pout = lastDiameter + 2*thisP(pIDX).d_o;
                                lastDiameter = d_Pout;
                                thisP(pIDX).diameter = (d_Pin + d_Pout)/2;
                                thisP(pIDX).length = pi*thisP(pIDX).diameter*thisP(pIDX).N;
                                
                                if isequal(thisP(pIDX).bifilar, 2)
                                    pIDX = pIDX + 2;
                                    thisP(pIDX + 1).diameter = thisP(pIDX).diameter;
                                    thisP(pIDX + 1).length = thisP(pIDX).length;
                                else
                                    pIDX = pIDX + 1;
                                end
                                
                                pCurrentLayer = 1;
                            end
                        else
                            if pCurrentLayer < thisP.N_L
                                if isequal(pCurrentLayer, 1)
                                    d_Pin = lastDiameter;
                                end
                                
                                lastDiameter = lastDiameter + 2*thisP.d_o;
                            else
                                d_Pout = lastDiameter + 2*thisP.d_o;
                                lastDiameter = d_Pout;
                                thisP.diameter = (d_Pin + d_Pout)/2;
                                thisP.length = pi*thisP.diameter*thisP.N;
                            end
                        end
                        
                        lastPS = 'P';
                    case 'P' % secondary layer
                        if nws > 1
                            if sCurrentLayer < thisS(sIDX).N_L
                                if isequal(sCurrentLayer, 1)
                                    d_Sin = lastDiameter;
                                end
                                
                                lastDiameter = lastDiameter + 2*thisS(sIDX).d_o;
                            else
                                d_Sout = lastDiameter + 2*thisS(sIDX).d_o;
                                lastDiameter = d_Sout;
                                thisS(sIDX).diameter = (d_Sin + d_Sout)/2;
                                thisS(sIDX).length = pi*thisS(sIDX).diameter*thisS(sIDX).N;
                                
                                if isequal(thisS(sIDX).bifilar, 2)
                                    sIDX = sIDX + 2;
                                    thisS(sIDX + 1).diameter = thisS(sIDX).diameter;
                                    thisS(sIDX + 1).length = thisS(sIDX).length;
                                else
                                    sIDX = sIDX + 1;
                                end
                                
                                sCurrentLayer = 1;
                            end
                        else
                            if sCurrentLayer < thisS.N_L
                                if isequal(sCurrentLayer, 1)
                                    d_Sin = lastDiameter;
                                end
                                
                                lastDiameter = lastDiameter + 2*thisS.d_o;
                            else
                                d_Sout = lastDiameter + 2*thisS.d_o;
                                lastDiameter = d_Sout;
                                thisP.diameter = (d_Pin + d_Pout)/2;
                                thisP.length = pi*thisP.diameter*thisP.N;
                            end
                        end
                        
                        lastPS = 'S';
                end
            end
            
            % quick check that winding height < window height
            if ~(lastDiameter < Core.d_center + 2*Core.window.height)
                warning('Windings exceed window height as arranged.')
            end
    end
end