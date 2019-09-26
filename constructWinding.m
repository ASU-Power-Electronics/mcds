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
    Wgap = Winding.t_ins;
    
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
                                  titleStr, 1, {'0', '0', '1'});
                
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
                              1, {'0', '0', '1'});

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
                                  titleStr, 1, {'0', '0', '1'});

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
                              1, {'0', '0', '1'});

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
    
    [thisP, thisS] = arrangeWindings(thisP, nwp, thisS, nws, Wgap, Core);

    Winding.primary = orderfields(thisP);
    Winding.secondary = orderfields(thisS);
    clear thisP thisS flag count bb isTrueMultiWinding isPaired
    clear isDirectPairing isNearDirectPair isOverallPairing isNearOverallPair
    clear isEligibleForInterleave

    % fun plot of core and winding geometry; basic validation of feasibility
    if strcmp(Core.c_leg_type, 'round')
        coreRadius = sqrt(Core.A_e/pi);
        figure
        polarplot(0:pi/180:2*pi, coreRadius*ones(1, 361), 'LineWidth', 2, 'Color', 'k');
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
        polarplot(0:pi/180:2*pi, outerRadius*ones(1, 361), 'LineWidth', 2, 'Color', 'k');
    else
        tBobbin = Core.window.height - Core.bobbin.height;
        hC = Core.d_center + tBobbin; % center leg height w/ bobbin (x)
        dC = Core.d_center2 + tBobbin; % center leg depth w/ bobbin (z)
        offsetX = hC/2;
        offsetY = dC/2;
        bifpair = false;
        
        figure
        rectangle('Position', [-offsetX, -offsetY, hC, dC], 'LineWidth', 2, 'EdgeColor', 'k')
        hold on
        
        if nwp > 1
            for idx = 1:nwp
                hW = Winding.primary(idx).N_L*Winding.primary(idx).d_o;
                rectangle('Position', [-offsetX - hW/2, -offsetY - hW/2, 2*offsetX + hW, 2*offsetY + hW], 'Curvature', 0.2)
                
                if isequal(Winding.primary(idx).bifilar, 2)
                    if bifpair
                        bifpair = false;
                        offsetX = offsetX + hW/2;
                        offsetY = offsetY + hW/2;
                    else
                        bifpair = true;
                    end
                else
                    offsetX = offsetX + hW/2;
                    offsetY = offsetY + hW/2;
                end
            end
        else
            hW = Winding.primary.N_L*Winding.primary.d_o;
            rectangle('Position', [-offsetX - hW/2, -offsetY - hW/2, 2*offsetX + hW, 2*offsetY + hW], 'Curvature', 0.2)
                
            if isequal(Winding.primary.bifilar, 2)
                bifpair = true;
            else
                offsetX = offsetX + hW/2;
                offsetY = offsetY + hW/2;
            end
        end
        
        if nws > 1
            for idx = 1:nws
                hW = Winding.secondary(idx).N_L*Winding.secondary(idx).d_o;
                rectangle('Position', [-offsetX - hW/2, -offsetY - hW/2, 2*offsetX + hW, 2*offsetY + hW], 'Curvature', 0.2)
                
                if isequal(Winding.primary(idx).bifilar, 2)
                    if bifpair
                        bifpair = false;
                        offsetX = offsetX + hW/2;
                        offsetY = offsetY + hW/2;
                    else
                        bifpair = true;
                    end
                else
                    offsetX = offsetX + hW/2;
                    offsetY = offsetY + hW/2;
                end
            end
        else
            hW = Winding.secondary.N_L*Winding.secondary.d_o;
            rectangle('Position', [-offsetX - hW/2, -offsetY - hW/2, 2*offsetX + hW, 2*offsetY + hW], 'Curvature', 0.2)
        end
        
        innerEdge = Core.window.height + Core.d_center/2;
        outerEdge = innerEdge + Core.d_center/2; % approximate
        rectangle('Position', [-outerEdge, -dC/2, outerEdge - innerEdge, dC], 'LineWidth', 2, 'EdgeColor', 'k')
        rectangle('Position', [innerEdge, -dC/2, outerEdge - innerEdge, dC], 'LineWidth', 2, 'EdgeColor', 'k')
    end
    
    hold off
    grid on
    axis equal

%     clear nwp nws coreRadius outerRadius
end