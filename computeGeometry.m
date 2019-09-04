%% computeGeometry
% Computes winding geometry from prompt response.

function Winding = computeGeometry(Winding, bb, answer)
    Winding.AWG_s = str2double(answer{1});
    Winding.NPW = str2double(answer{3});
    Winding.N_s = str2double(answer{2})*Winding.NPW;

    [Winding.d_s, Winding.A_s] = AWG2m(Winding.AWG_s);
    Winding.d_o = sqrt(Winding.N_s/Winding.NPW)*Winding.d_s/0.68125;
    Winding.tplMax = floor(bb/(Winding.NPW*Winding.bifilar*Winding.d_o));
    [Winding.tpl, Winding.N_L] = conjureLayers(Winding);
    Winding.A_w_Cu = Winding.A_s*Winding.N_s;
    Winding.A_Cu = Winding.A_w_Cu*Winding.N;
    Winding.AWG_e = m2AWG(Winding.A_w_Cu);
end