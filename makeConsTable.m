%% makeConsTable
% Creates a table of constructions given an existing (or empty) table.

function T = makeConsTable(T, cons)
    nRows = size(cons, 2);
    varNames = {'Bundle_AWG', 'Bundles_per_Layer'};
    rowNames = cell(1, nRows);

    for row = 1:nRows
        name = sprintf('%d/%d', cons(row).N_s, cons(row).AWG);
        rowNames{1, row} = name;
    end
       
    if isempty(T)        
        T = table([cons(:).AWG_e]', [cons(:).N_bpl]', ...
                  'VariableNames', varNames, ...
                  'RowNames', rowNames);
    else
        Tadd = table([cons(:).AWG_e]', [cons(:).N_bpl]', ...
                     'VariableNames', varNames, ...
                     'RowNames', rowNames);
        T = [T; Tadd];
    end
end