%% constructWinding
% Constructs a transformer winding, calculating layer geometry using selected
% wires.

function Winding = constructWinding(Transformer)
    Winding = Transformer.winding;
    Properties = Transformer.properties;
    Core = Transformer.core;
    thisP = Winding.primary;
    thisS = Winding.secondary;
    nwp = Properties.N_wp;
    nws = Properties.N_ws;
    bb = Core.bobbin.breadth;
    Winding.A_Cu = 0;
    
    % add fields to all primary windings
    Z = num2cell(zeros(size(thisP(:))));
    [thisP(:).bifilar] = Z{:};
    [thisP(:).AWG_s] = Z{:};
    [thisP(:).NPW] = Z{:};
    [thisP(:).N_s] = Z{:};
    [thisP(:).d_s] = Z{:};
    [thisP(:).A_s] = Z{:};
    [thisP(:).d_o] = Z{:};
    [thisP(:).tplMax] = Z{:};
    [thisP(:).tpl] = Z{:};
    [thisP(:).N_L] = Z{:};
    [thisP(:).A_w_Cu] = Z{:};
    [thisP(:).A_Cu] = Z{:};
    [thisP(:).AWG_e] = Z{:};
    [thisP(:).diameter] = Z{:};
    [thisP(:).length] = Z{:};
    
    % add fields to all secondary windings
    Z = num2cell(zeros(size(thisS(:))));
    [thisS(:).bifilar] = Z{:};
    [thisS(:).AWG_s] = Z{:};
    [thisS(:).NPW] = Z{:};
    [thisS(:).N_s] = Z{:};
    [thisS(:).d_s] = Z{:};
    [thisS(:).A_s] = Z{:};
    [thisS(:).d_o] = Z{:};
    [thisS(:).tplMax] = Z{:};
    [thisS(:).tpl] = Z{:};
    [thisS(:).N_L] = Z{:};
    [thisS(:).A_w_Cu] = Z{:};
    [thisS(:).A_Cu] = Z{:};
    [thisS(:).AWG_e] = Z{:};
    [thisS(:).diameter] = Z{:};
    [thisS(:).length] = Z{:};
    
    fprintf('\nWire Selection:\n')

    % primary winding(s) layer geometry
    if nwp > 1
        for idx = 1:nwp
            bifAnswer = '';

            while isequal(bifAnswer, '')
                titleStr = ['Primary Winding ', num2str(idx), ' Wire Selection'];
                answer = inputdlg({'Wire gauge', ...
                                   'Number of strands', ...
                                   'Number of parallel wires'}, ...
                                  titleStr, 1, {'0', '0', '0'});
                
                bifAnswer = questdlg('Bifilar-wound?');

                if isequal(bifAnswer, 'Yes')
                    thisP(idx).bifilar = 2;
                elseif isequal(bifAnswer, 'No')
                    thisP(idx).bifilar = 1;
                elseif isequal(bifAnswer, 'Cancel') || isequal(bifAnswer, '')
                    bifAnswer = '';
                end
            end
            
            thisP(idx) = computeGeometry(thisP(idx), bb, answer);

            Winding.A_Cu = Winding.A_Cu + thisP(idx).A_Cu;

            fprintf('Primary %d strand AWG:  %d\n', idx, thisP(idx).AWG_s)
            fprintf('Primary %d strand count: %d\n', idx, thisP(idx).N_s)
            fprintf('Primary %d equivalent AWG: %d\n', idx, round(thisP(idx).AWG_e))
            fprintf('Primary %d turns per layer:  %d\n', idx, thisP(idx).tpl)
            fprintf('Primary %d layer count:  %d\n', idx, thisP(idx).N_L)
        end
    else
        bifAnswer = '';

        while isequal(bifAnswer, '')
            answer = inputdlg({'Wire gauge', ...
                               'Number of strands', ...
                               'Number of parallel wires'}, ...
                              'Primary Winding Wire Selection', ...
                              1, {'0', '0', '0'});

            bifAnswer = questdlg('Bifilar-wound?');

            if isequal(bifAnswer, 'Yes')
                thisP.bifilar = 2;
            elseif isequal(bifAnswer, 'No')
                thisP.bifilar = 1;
            elseif isequal(bifAnswer, 'Cancel') || isequal(bifAnswer, '')
                bifAnswer = '';
            end
        end

        thisP = computeGeometry(thisP, bb, answer);

        Winding.A_Cu = Winding.A_Cu + thisP.A_Cu;

        fprintf('Primary strand AWG:  %d\n', thisP.AWG_s)
        fprintf('Primary strand count: %d\n', thisP.N_s)
        fprintf('Primary equivalent AWG: %d\n', round(thisP.AWG_e))
        fprintf('Primary turns per layer:  %d\n', thisP.tpl)
        fprintf('Primary layer count:  %d\n', thisP.N_L)
    end

    % secondary winding(s) layer geometry
    if nws > 1
        for idx = 1:nws
            bifAnswer = '';

            while isequal(bifAnswer, '')
                titleStr = ['Secondary Winding ', num2str(idx), ' Wire Selection'];
                answer = inputdlg({'Wire gauge', ...
                                   'Number of strands', ...
                                   'Number of parallel wires'}, ...
                                  titleStr, 1, {'0', '0', '0'});

                bifAnswer = questdlg('Bifilar-wound?');

                if isequal(bifAnswer, 'Yes')
                    thisS(idx).bifilar = 2;
                elseif isequal(bifAnswer, 'No')
                    thisS(idx).bifilar = 1;
                elseif isequal(bifAnswer, 'Cancel') || isequal(bifAnswer, '')
                    bifAnswer = '';
                end
            end

            thisS(idx) = computeGeometry(thisS(idx), bb, answer);

            Winding.A_Cu = Winding.A_Cu + thisS(idx).A_Cu;

            fprintf('Secondary %d strand AWG:  %d\n', idx, thisS(idx).AWG_s)
            fprintf('Secondary %d strand count: %d\n', idx, thisS(idx).N_s)
            fprintf('Secondary %d equivalent AWG: %d\n', idx, round(thisS(idx).AWG_e))
            fprintf('Secondary %d turns per layer:  %d\n', idx, thisS(idx).tpl)
            fprintf('Secondary %d layer count:  %d\n', idx, thisS(idx).N_L)
        end
    else
        bifAnswer = '';

        while isequal(bifAnswer, '')
            answer = inputdlg({'Wire gauge', ...
                               'Number of strands', ...
                               'Number of parallel wires'}, ...
                              'Secondary Winding Wire Selection', ...
                              1, {'0', '0', '0'});

            bifAnswer = questdlg('Bifilar-wound?');

            if isequal(bifAnswer, 'Yes')
                thisS.bifilar = 2;
            elseif isequal(bifAnswer, 'No')
                thisS.bifilar = 1;
            elseif isequal(bifAnswer, 'Cancel') || isequal(bifAnswer, '')
                bifAnswer = '';
            end
        end

        thisS = computeGeometry(thisS, bb, answer);

        Winding.A_Cu = Winding.A_Cu + thisS.A_Cu;

        fprintf('Secondary strand AWG:  %d\n', thisS.AWG_s)
        fprintf('Secondary strand count: %d\n', thisS.N_s)
        fprintf('Secondary equivalent AWG: %d\n', round(thisS.AWG_e))
        fprintf('Secondary turns per layer:  %d\n', thisS.tpl)
        fprintf('Secondary layer count:  %d\n', thisS.N_L)
    end
    
    [thisP, thisS] = arrangeWindings(thisP, nwp, thisS, nws, Core);

    Winding.primary = orderfields(thisP);
    Winding.secondary = orderfields(thisS);
    clear thisP thisS flag count bb isTrueMultiWinding isPaired
    clear isDirectPairing isNearDirectPair isOverallPairing isNearOverallPair
    clear isEligibleForInterleave

    % fun plot of core and winding geometry; basic validation of feasibility,
    % though for rectangular center legs, it can get a bit wonky
    coreRadius = sqrt(Core.A_e/pi);
    figure;
    polarplot(0:pi/180:2*pi, coreRadius*ones(1, 361));
    hold on;

    if nwp > 1
        for idx = 1:nwp
            r = Winding.primary(idx).diameter/2;
            polarplot(0:pi/180:2*pi, r*ones(1, 361));
        end
    else
        r = Winding.primary.diameter/2;
        polarplot(0:pi/180:2*pi, r*ones(1, 361));
    end

    if nws > 1
        for idx = 1:nws
            r = Winding.secondary(idx).diameter/2;
            polarplot(0:pi/180:2*pi, r*ones(1, 361));
        end
    else
        r = Winding.secondary.diameter/2;
        polarplot(0:pi/180:2*pi, r*ones(1, 361));
    end

    outerRadius = Core.window.height + Core.d_center/2;
    polarplot(0:pi/180:2*pi, outerRadius*ones(1, 361));
    hold off
    grid on

    clear nwp nws coreRadius outerRadius
end