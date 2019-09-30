%% computeTempRise
% Computes temperature rise of transformer based on very simple thermal model
% given by Kazimierczuk, using computed loss quantities and datasheet properties
% of core, alongside material coefficients.
%
% kc is given as 3.5 to 5 for both MnZn and NiZn in Ferroxcube Data Handbook
% kw is given as 1/2.51 for a particular litz wire in Biela; the value here
% is just a temporary value for testing
% thetaSA is from Hurley (1998), given as ~1/10 with backing from 2 sources

function [T, DT] = computeTempRise(P, C, W)
    TA = 25; % Ambient temperature [°C]
    
    % Core
    rc = sqrt(C.A_e/pi); % Effective core radius [m]
    kc = 5; % Core thermal conductivity [W/(m°C)]
    Asc = C.A_s; % Effective core surface area [m^2]
    thetaC = rc/(kc*Asc); % Core thermal resistance [°C/W]
    
    % Primary
    twp = [W.primary.N_L]*[W.primary.d_o]'; % Primary winding thickness (height) [m]
    kwp = 1/2.5; % Primary winding thermal conductivity [W/(m°C)]
    Aswp = pi*([W.primary.NPW] + [W.primary.tpl] - 1)*[W.primary.d_o]'; % Primary winding surface area [m^2]
    thetaP = twp/(kwp*Aswp); % Primary winding(s) thermal resistance [°C/W]
    
    % Secondary
    tws = [W.secondary.N_L]*[W.secondary.d_o]'; % Secondary winding(s) thickness (height) [m]
    kws = 1/2.5; % Secondary winding(s) thermal conductivity [W/(m°C)]
    Asws = pi*([W.secondary.NPW] + [W.secondary.tpl] - 1)*[W.secondary.d_o]'; % Secondary winding(s) surface area [m^2]
    thetaS = tws/(kws*Asws); % Secondary winding(s) thermal resistance [°C/W]
    thetaSA = 1/10; % Secondary winding(s) to ambient thermal resistance [°C/W]
    
    % Loss values
    Pfe = P.P_Fe; % Core loss [W]
    PCuP = sum([W.primary.P_Cu]); % Primary winding loss [W]
    PCuS = sum([W.secondary.P_Cu]); % Secondary winding loss [W]
    
    % Final Temperatures
    DT = thetaC*Pfe + thetaP*(Pfe + PCuP) + (thetaS + thetaSA)*(Pfe + PCuP + PCuS); % Temperature rise [°C]
    T = TA + DT; % Final surface temperature [°C]
end