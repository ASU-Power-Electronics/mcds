%% coreMaterial
% Selects power-oriented ferrite for transformer design by switching frequency;
% offers option for manual material selection by user input.

function material = coreMaterial(f, Bpk)
    % Table of materials used (Ferroxcube)
    % NOTE: low-temp B_sat in this table, 100 deg C used in file
    % material	mu_i	B_sat	T_c	rho	ferriteType	f_1	f_2  f_c
    % ~~~~3C85      2000    400     200 2   MnZn        0   2e5  1e5~~ removed
    % 3C90      2300	470     220	5	MnZn        0	2e5  1e5
    % 3C92      1500	520     280	5	MnZn        0	2e5  1e5
    % 3C91      3000	470     220	5	MnZn        0	3e5  15e4
    % 3C94      2300	470     220	5	MnZn        0	3e5  15e4
    % 3C96      2000	500     240	5	MnZn        0	4e5  2e5
    % ~~~~3C98      2500	530     230	5	MnZn        0	4e5  2e5~~ removed
    % 3C93      1800	500     240	5	MnZn        0	5e5  25e4
    % 3C95      3000	530     215	5	MnZn        0	5e5  25e4
    % 3C97      3000	530     215	5	MnZn        0	5e5  25e4
    % 3F3       2000	440     200	2	MnZn        2e5	5e5  35e4
    % 3F31      1800	520     230	8	MnZn        2e5	5e5  35e4
    % 3F35      1400	500     240	10	MnZn        5e5	1e6  75e4
    % 3F36      1600	520     230	12	MnZn        5e5	1e6  75e4
    % 3F4       900     410     220	10	MnZn        1e6	2e6  15e5
    % 3F45      900     420     300	10	MnZn        1e6	2e6  15e5
    % 3F5       650     380     300	10	MnZn        2e6	4e6  3e6
    % 4F1       80      320     260	1e5	NiZn        4e6	15e6 95e5

    ok = 0;
    answer = {};

    while isequal(ok, 0)
        mode = questdlg('Automatic or manual material selection?', ...
                        'Material selection', ...
                        'Automatic', ...
                        'Manual', ...
                        'Automatic');

        load('Materials.mat', 'Materials');
        [~, numMats] = size(Materials);
        fVec = [Materials.f_c]; % populate center frequencies
        f2 = [Materials.f_2]; % populate upper bound
        fMax = max(f2);
        material = struct;

        switch mode
            case 'Automatic'
                if f < fMax % ensure material exists for frequency
                    fDist = abs(fVec - f); % determine 1-norm distances to center
                    fSpread = abs(fVec - f2); % determine 1-norm spreads
                    c = find(fDist < fSpread); % determine indices of fits
                    
                    % reduce set by saturation flux density (~85%, since value
                    % is approximate)
                    c([Materials(c).B_sat] < 0.85*Bpk) = [];

                    muVec = zeros(1, length(c));
                    
                    for opt = 1:length(c) % iterate through options to populate
                        material.opts(opt) = Materials(c(opt));
                        mu_r = evaluatePermeability(material.opts(opt).mu_a, Bpk);
                        muVec(opt) = mu_r;
                    end

                    muMax = max(muVec);
                    idx = find(muVec == muMax);
                    
                    if length(idx) == 1 % if only one choice
                        material.main = material.opts(idx);
                    else % if multiple choices, pick first
                        material.main = material.opts(idx(1));
                    end
                    
                    ok = 1;

                else % no can do, buckaroo
                    warning('No core material available for the selected frequency.')
                end
            case 'Manual'
                while isequal(answer, {})
                    [Selection, ok] = listdlg('ListString', ...
                                              {'Add New', ...
                                               '3C90', ...
                                               '3C92', ...
                                               '3C91', ...
                                               '3C94', ...
                                               '3C96', ...
                                               '3C93', ...
                                               '3C95', ...
                                               '3C97', ...
                                               '3F3', ...
                                               '3F31', ...
                                               '3F35', ...
                                               '3F36', ...
                                               '3F4', ...
                                               '3F45', ...
                                               '3F5', ...
                                               '4F1'}, ...
                                              'SelectionMode', 'single', ...
                                              'Name', 'Ferrite Material Selection');

                    if isequal(Selection, 1)
                        opts = struct;
                        opts.Interpreter = 'tex';
                        opts.Resize = 'on';
                        answer = inputdlg({'Material name', ...
                                           'Initial Permeability (\mu_i)', ...
                                           'Amplitude Permeability (\mu_a)', ...
                                           'Saturation Flux Density (B_sat)', ...
                                           'Curie Temperature (T_c)', ...
                                           'Resistivity (\rho)', ...
                                           'Ferrite Type (ex. MnZn, NiZn)', ...
                                           'Frequency Lower Bound (in Hz)', ...
                                           'Frequency Upper Bound'}, ...
                                          'Add new material', 1, ...
                                          {'0', '0', '0', '0', 'None', '0', num2str(fMax)}, opts);
                        material.opts(1).material = answer{1};
                        material.opts(1).mu_i = str2double(answer{2});
                        material.opts(1).mu_a = str2num(answer{3}); % hopefully non-scalar
                        material.opts(1).B_sat = str2double(answer{4});
                        material.opts(1).T_c = str2double(answer{5});
                        material.opts(1).rho = str2double(answer{6});
                        material.opts(1).ferriteType = answer{7};
                        material.opts(1).f_1 = str2double(answer{8});
                        material.opts(1).f_2 = str2double(answer{9});
                        material.opts(1).f_c = (material.opts(1).f_1 + material(1).opts.f_2)/2;
                        material.main = material.opts(1);
                        
                        % write back to mat file (not actually looping...)
                        Materials(numMats + 1) = material.main;
                        save('Materials.mat', 'Materials')
                    elseif ~isequal(ok, 0) % selection made
                        material.opts = Materials(Selection - 1);
                        material.main = material.opts;
                        answer = {1};
                    else % canceled manual input
                        answer = {1};
                    end
                end
        end
    end
    
    clear Materials
end