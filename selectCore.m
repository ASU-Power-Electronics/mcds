%% selectCore
% Selects a ferrite core for use with transformer design, given base core
% parameters, base winding parameters, and base transformer properties.  Core
% geometry coefficient values in the associated database assume a window
% utilization factor of 0.3.  MLT values are calculated based on mid-window
% diamater for round center legs, and window height + center leg distance for
% rectangular cores in both x and y directions, for 4 segments equalling one
% full turn passing through the center of the window.
%TODO: find existing core if name/params are the same and edit

function Core = selectCore(Core, Properties)
    material = Core.material; % placeholder
    matName = material.main.material;
    coreNames = {};
    load('Cores.mat', 'Cores')
    
    % cull cores for which material is not available
    [~, numCores] = size(Cores);
    indices = [];
    count = 0;
    
    for core = 1:numCores
        if any(strcmp(Cores(core).materials, matName))
            count = count + 1;
            indices = [indices, core];
            coreNames{count} = Cores(core).name;
        end
    end
    
    CoresShort = Cores(indices);
    
    numCoresOrig = numCores;
    [~, numCores] = size(CoresShort);
    KgList = zeros(1, numCores);
    
    % clean up structure
    CoresShort = rmfield(CoresShort, 'materials');
    
    % extract Kg values
    for core = 1:numCores
        KgList(core) = CoresShort(core).K_g;
    end
    
    % prepare minimum value for comparison
    KgMin = Properties.K_g;
    
    answer = questdlg('Automatic or manual core selection?', ...
                      'Core Selection', ...
                      'Automatic', ...
                      'Manual', ...
                      'Automatic');
    
    switch answer
        case 'Automatic' % find tightest Kg fit, select options within 10% Kg
            idx = [];
            KgFit = KgList - KgMin;
            minKgFit = min(KgFit(KgFit > 0));
            
            if isempty(minKgFit)
                minKgFit = -min(abs(KgFit));
                warning('No cores in database with Kg greater than minimum.  Results may be suboptimal.')
            end
            
            [~, best] = find(KgFit == minKgFit);
            Core = CoresShort(best);
            Core.material = material;
            
            for core = best:numCores
                if CoresShort(core).K_g <= 1.1*Core.K_g
                    idx = [idx, core];
                end
            end
            
            if length(idx) > 1
                for core = 1:length(idx)
                    Core.opts(core) = CoresShort(idx(core));
                end
            else
                Core.opts = orderfields(Core);
            end
        case 'Manual' % select from list
            ok = 0;
            coreNames = [{'Add New'}, coreNames];
            
            while isequal(ok, 0)
                [Selection, ok] = listdlg('ListString', ...
                                          coreNames, ...
                                          'SelectionMode', 'single', ...
                                          'Name', 'Core Selection');
            
                if ok && ~isequal(Selection, 1)
                	Core = CoresShort(Selection - 1);
                    
                    if Core.K_g < KgMin
                        resp = questdlg('Core does not meet minimum core geometry coefficient, choose a larger core?', 'Core Geometry Coefficient Warning');
                        
                        if strcmp(resp, 'Yes')
                            ok = 0;
                            continue
                        elseif strcmp(resp, 'No')
                            ok = 1;
                        else
                            warning("Core does not meet minimum core geometry coefficient, choose a larger core.")
                            ok = 0;
                            break
                        end
                    else
                        ok = 1;
                    end
                    
                    Core.opts = CoresShort([1:Selection - 2, Selection:end]);
                    Core.material = material;
                elseif ok && isequal(Selection, 1)
                    coreData = inputdlg({'Core designation (ex. ER32_6_25)', ...
                                         'Center leg type (round, rectangular)', ...
                                         'Effective volume [mm^3]', ...
                                         'Effective length [mm]', ...
                                         'Effective area [mm^2]', ...
                                         'Bobbin breadth [mm] (optional)', ...
                                         'Bobbin diameter [mm] (optional)'}, ...
                                        'Core Data Input', 1, ...
                                        {'none', 'round/rectangular', '0', '0', '0', '', ''});
                    [geometry, As, measures] = getGeometry(coreData{1});
                    
                    % store data
                    Core.name = coreData{1};
                    Core.c_leg_type = coreData{2};
                    Core.V_e = str2double(coreData{3})*1e-9;
                    Core.l_e = str2double(coreData{4})*1e-3;
                    Core.A_e = str2double(coreData{5})*1e-6;
                    Core.A_s = As;
                    Core.geometry = geometry;
                    
                    % extract and/or compute
                    Core.d_center = measures.d_center;
                    Core.d_center2 = measures.d_center2;
                    Core.window = measures.window;

                    % calculate remaining data
                    if ~isequal(coreData{10}, '')
                        Core.bobbin.breadth = str2double(coreData{10})*1e-3;
                    else
                        Core.bobbin.breadth = Core.window.breadth - 1e-3;
                    end

                    if ~isequal(coreData{11}, '')
                        Core.bobbin.height = str2double(coreData{11})*1e-3;
                    else
                        Core.bobbin.height = Core.window.height - 0.5e-3;
                    end

                    Core.W_a = Core.bobbin.breadth*Core.bobbin.height;

                    if strcmp(coreData{2}, 'round')
                        Core.MLT = pi*(Core.d_center + Core.window.height);
                    else
                        offset = Core.window.height/2;
                        Core.MLT = 2*(Core.d_center + Core.d_center2) + 8*offset;
                    end

                    Core.A_p = Core.A_e*Core.W_a;
                    Core.K_g = Core.A_p*Core.A_e*0.3/Core.MLT; % K_u = 0.3
                    
                    if Core.K_g < KgMin
                        resp = questdlg('Core does not meet minimum core geometry coefficient, choose a larger core?', 'Core Geometry Coefficient Warning');
                        
                        if strcmp(resp, 'Yes')
                            ok = 0;
                            continue
                        elseif strcmp(resp, 'No')
                            ok = 1;
                        else
                            warning("Core does not meet minimum core geometry coefficient, choose a larger core.")
                            ok = 0;
                            break
                        end
                    else
                        ok = 1;
                    end
                    
                    coreSaver = Core;
                    Core = orderfields(Core);
                    Core.opts = Core;
                    
                    % save to structure after swapping material defs
                    saveIt = questdlg('Save this core?', 'Save Core');
                    
                    if isequal(saveIt, 'Yes')
                        % check for existing core
                        if ismember(coreData{1}, [Cores(:).name])
                            idx = strcmp({Cores(:).name}, coreData{1});
                            
                            if ~ismember(Core.material.main.material, Cores(idx).materials(:))
                                Cores(idx).materials = [Cores(idx).materials(:)', {Core.material.main.material}];
                            end
                        else
                            coreSaver.materials = {Core.material.main.material};
                            coreSaver = rmfield(coreSaver, 'material');
                            coreSaver.c_leg_type = coreData{2};
                            Cores(numCoresOrig + 1) = coreSaver;
                        end
                        
                        save('Cores.mat', 'Cores')
                    end
                end
            end
            
            % Offer to start over if cancelled, return results as nested call
            % Otherwise crash (aka "exit")
            if isequal(ok, 0)
                resp = questdlg('Select New Core or Exit?', ...
                                'Continue or Exit', ...
                                'New Core', 'Exit', 'New Core');
                            
                if strcmp(resp, 'New Core')
                    Core = selectCore(Core, Properties);
                end                    
            end
    end
    
    clear Cores
end