%% selectCore
% Selects a ferrite core for use with transformer design, given base core
% parameters, base winding parameters, and base transformer properties.  Core
% geometry coefficient values in the associated database assume a window
% utilization factor of 0.3.  MLT values are calculated based on mid-window
% diamater for round center legs, and window height + center leg distance for
% rectangular cores in both x and y directions, for 4 segments equalling one
% full turn passing through the center of the window.

%TODO: change from warning to message box

function Core = selectCore(Core, Properties)
    material = Core.material; % placeholder
    matName = material.main.material;
    coreNames = {};
    load('Cores.mat', 'Cores')
    
    % cull cores for which material is not available
    [~, numCores] = size(Cores);
    indices = [];
    count = 0;
    
    for i = 1:numCores
        if any(strcmp(Cores(i).materials, matName))
            count = count + 1;
            indices = [indices, i];
            coreNames{count} = Cores(i).name;
        end
    end
    
    CoresShort = Cores(indices);
    
    numCoresOrig = numCores;
    [~, numCores] = size(CoresShort);
    KgList = zeros(1, numCores);
    
    % clean up structure
    CoresShort = rmfield(CoresShort, 'materials');
    
    % extract Kg values
    for i = 1:numCores
        KgList(i) = CoresShort(i).K_g;
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
            [~, best] = find(KgFit == minKgFit);
            Core = CoresShort(best);
            Core.material = material;
            
            for i = best:numCores
                if CoresShort(i).K_g <= 1.1*Core.K_g
                    idx = [idx, i];
                end
            end
            
            if length(idx) > 1
                for i = 1:length(idx)
                    Core.opts(i) = CoresShort(idx(i));
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
                        warning("Core does not meet minimum core geometry coefficient, choose a larger core.")
                        ok = 0;
                        continue
                    else
                        ok = 1;
                    end
                    
                    Core.opts = CoresShort([1:Selection - 2, Selection:end]);
                    Core.material = material;
                elseif ok && isequal(Selection, 1)
                    coreSaver = struct();
                    coreData = inputdlg({'Core designation (ex. ER32_6_25)', ...
                                         'Center leg type (round, rectangular)', ...
                                         'Effective volume [mm^3]', ...
                                         'Effective length [mm]', ...
                                         'Effective area [mm^2]', ...
                                         'Center leg height [mm]', ...
                                         'Center leg depth [mm]', ...
                                         'Window breadth [mm]', ...
                                         'Window height [mm]', ...
                                         'Bobbin breadth [mm] (optional)', ...
                                         'Bobbin height [mm] (optional)'}, ...
                                        'Core Data Input', 1, ...
                                        {'none', '0', '0', '0', '0', '0', '0', '0', '0', '', ''});
                    % store data
                    Core.name = coreData{1};
                    Core.V_e = str2double(coreData{3})*1e-9;
                    Core.l_e = str2double(coreData{4})*1e-3;
                    Core.A_e = str2double(coreData{5})*1e-6;
                    Core.d_center = str2double(coreData{6})*1e-3;
                    Core.d_center2 = str2double(coreData{7})*1e-3;
                    Core.window.breadth = str2double(coreData{8})*1e-3;
                    Core.window.height = str2double(coreData{9})*1e-3;
                    
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
                    
                    if isequal(coreData{2}, 'round')
                        Core.MLT = pi*(Core.d_center + Core.window.height);
                    else
                        offset = Core.window.height/2;
                        Core.MLT = 2*(Core.d_center + Core.d_center2) + 8*offset;
                    end
                    
                    Core.A_p = Core.A_e*Core.W_a;
                    Core.K_g = Core.A_p*Core.A_e*0.3/Core.MLT; % K_u = 0.3
                    
                    if Core.K_g < KgMin
                        warning("Core does not meet minimum core geometry coefficient, choose a larger core.")
                        ok = 0;
                        continue
                    else
                        ok = 1;
                    end
                    
                    coreSaver = Core;
                    Core = orderfields(Core);
                    Core.opts = Core;
                    
                    % save to structure after swapping material defs
                    saveIt = questdlg('Save this core?', 'Save Core');
                    
                    if isequal(saveIt, 'Yes')
                        coreSaver.materials = {Core.material.main.material};
                        coreSaver = rmfield(coreSaver, 'material');
                        coreSaver.c_leg_type = coreData{2};
                        Cores(numCoresOrig + 1) = coreSaver;
                        save('Cores.mat', 'Cores')
                    end
                end
            end
    end
    
    clear Cores
end