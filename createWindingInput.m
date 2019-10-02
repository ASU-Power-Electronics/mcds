%% createWindingInput
% Creates an input dialog for a winding structure, cardinality-aware.  Input
% structure has fields:
%
% prompt: 2 x 1 or 2 x 2 depending on where winding number is placed (end or
%         middle, respectively)
% dlgtitle: character vector
% definput: scalar or 2 x 1 vector of default values
% np: number of primary windings
% ns: number of secondary windings

function response = createWindingInput(inputStruct)
    promptBase = inputStruct.prompt;
    dlgtitle = inputStruct.dlgtitle;
    definputBase = inputStruct.definput;
    np = inputStruct.np;
    ns = inputStruct.ns;
    prompt = {};
    definput = {};
    
    for p = 1:np
        if np > 1
            numstr = [' ', num2str(p), ' '];
        else
            numstr = ' ';
        end
        
        if isequal(size(promptBase, 2), 2)            
            prompt{1, end + 1} = [promptBase{1, 1}, numstr, promptBase{1, 2}];
        else
            prompt{1, end + 1} = [promptBase{1, 1}, numstr];
        end
        
        if isequal(size(definputBase, 1), 2)
            definput{1, end + 1} = num2str(definputBase{1, 1});
        else
            definput{1, end + 1} = num2str(definputBase);
        end
    end
    
    for s = 1:ns
        if ns > 1
            numstr = [' ', num2str(s), ' '];
        else
            numstr = ' ';
        end
        
        if isequal(size(promptBase, 2), 2)            
            prompt{1, end + 1} = [promptBase{2, 1}, numstr, promptBase{2, 2}];
        else
            prompt{1, end + 1} = [promptBase{2, 1}, numstr];
        end
        
        if isequal(size(definputBase, 1), 2)
            definput{1, end + 1} = num2str(definputBase{2, 1});
        else
            definput{1, end + 1} = num2str(definputBase);
        end
    end
    
    dims = 1;
    response = inputdlg(prompt, dlgtitle, dims, definput, 'on');
end