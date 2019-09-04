%% conjureLayers
% Conjures up a layer structure, maximizing number of turns per layer first.
% Returns the turns per layer and number of layers.

function [tpl, N_L] = conjureLayers(Winding)
    flag = 0;
    layers = 1;
    turns = Winding.N;

    while ~flag
        if isequal(mod(turns, layers), 0) && turns/layers <= Winding.tplMax
            tpl = turns/layers;
            N_L = layers;
            flag = 1;
        else
            layers = layers + 1;
        end
    end    
end