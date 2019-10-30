%% Transformer Design Script
% *R. Scott Mongrain*
% - _June 2016-October 2019_
%
% Transferred to ASU Power Electronics organization October 2019.
%
% This program is a tool for use in designing power electronics magnetic
% components (specifically transformers of various descriptions).
%
% In round figures, the design process goes as follows:
% 
% # Design converter with ideal transformer (PLECS preferred)
% # Provide design values and waveforms
% # Select core manually or automatically with guidance
% # Obtain LitzOpt input and run
% # Select winding configuration with guidance
% # Obtain necessary modifications to Ap/Kg
% # Iterate over the previous seven steps until satisfied with delta
% # Simulate (circuit and FEM), build, test, enjoy

%% Development Notes
% * Eventually would like to connect this to converterDesign.m for a more
% complete solution to that problem.
% * Inductor design needs a solution
% * GUI design?
% * *Transformer Structure Implementation*
% * - Can combine thisC, thisW, etc. into this as substructures (e.g. this.C)
% * *Other Stuff/Current Work*
% * - Implement iteration and pass counter token
% * - Implement D_eff for SFDT
% * - Generate Maxwell vbscripts based on core/windings
% * - Attach to PLECS via simulink blockset for faster iterations
% * - Investigate *webread* for scraping HTML
% * - https://blogs.mathworks.com/loren/2017/07/10/web-scraping-and-mining-unstructured-data-with-matlab/
% * - Better file system interface
% * - TODO's
% * *Testing*
% * - Compare number of iterations, tic/toc, etc. for this, direct A_p, and
% direct K_g
% * - Unit tests for associated modules including edge cases
% * - Test and validation of actual transformer designs

%% Program Setup
% This section clears all variables (but leaves other workspace elements alone),
% closes any open handles, cleans out the command window, and sets up the
% display of variables to the console.  It also establishes a tolerance to be
% used in any numerical methods employed.  Finally, it adds function libraries
% to the path necessary to the computation of several math functions.

clear variables
close all
clc
format shortEng
format compact

tol = eps*1e5; % ~2e-11 for 2-digit convergence at nano-scale

% https://www.mathworks.com/matlabcentral/fileexchange/37302-struve-functions
addpath('Struve01') % for inductance calculations

%% Constants and Structures Initializations
% This section initializes physical constants and device structures to be used
% throughout the problem.  As such, the variables are defined globally.
%%
% *Physical Constants*
%
% * $\mu_0$:  vacuum permeability in [H/m]
% * $\mu_{r, Cu}$:  copper relative permeability
% * $\mu_{Cu}$:  copper permeability in [H/m]
% * $\rho_{Cu}$:  copper resistivity in [Ohm*m]
% * $\sigma_{Cu}$:  copper conductivity in [S/m]
% * $J_{Max}$:  maximum current density in [A/m^2]
%
% NOTE:  All copper is treated as annealed copper.

global MU_0 MU_R_CU MU_CU RHO_CU SIGMA_CU J_MAX

MU_0 = 4*pi*1e-7;
MU_R_CU = 0.999994;
MU_CU = MU_R_CU*MU_0;
RHO_CU = 1.724e-8;
SIGMA_CU = 1/RHO_CU;
J_MAX = 6.8e6; % From fit to data in Odendaal/Ferreira (1999)

%%
% *Structure Initializations*
% 
% Time contains all time vectors used, at each time resolution.
% Converter contains information about the converter for which the transformer
% was designed.
% Transformer contains three substructures:
% - core contains information about the core, including material, geometry, and
%   some waveforms.
% - winding contains information about the windings, including construction and
%   geometry information, as well as all of the primary waveforms of converter
%   operation.
% - properties is a catch-all for transformer properties that do not
%   specifically pertain to only the core or windings.

Time = struct;
Converter = struct;
Transformer = struct('core', struct,...
                     'winding', struct,...
                     'properties', struct);

%% User Input and Design Specifications
% 
% This section accepts the user input (e.g. converter-specific values,
% transformer design specs, etc.), and other, temporary variables to use as
% initial points.  These values establish the physical parameters of the
% problem, including all design-specific variables related to the converter and
% transformer.  Now allows for an input frequency component for AC input.
%%
% *Converter Base Parameters*
% 
% * $V_{in, min}$:  minimum converter input voltage in [V]
% * $V_{in, max}$:  maximum converter input voltage in [V]
% * $V_o$:  nominal converter output voltage in [V]
% * $P_o$:  nominal converter output power in [W]
% * $f_s$:  converter switching frequency in [Hz]
%TODO: better file system interface starts here

button = questdlg('Create new converter or use existing?', ...
                  'Converter Selection', ...
                  'New', ...
                  'Existing', ...
                  'Existing');
switch button
    case 'New' % creates new converter from user input
        answer = inputdlg({'Name', ...
                           'Minimum Input Voltage [V]', ...
                           'Maximum Input Voltage [V]', ...
                           'Nominal Output Voltage [V]', ...
                           'Nominal Output Power [W]', ...
                           'Input AC Frequency (0 if DC) [Hz]', ...
                           'Switching Frequency [Hz]', ...
                           'Duty Ratio [#]'}, ...
                          'Converter Parameters', ...
                          1, ...
                          {'0', '0', '0', '0', '0', '0', '0', '0'}, 'on');
                      
        if isempty(answer) % cancel button handling
            error('Aborted, exiting design script.')
        end
        
        this = Converter; % one-way (get) alias for brevity
        this.name = answer{1};
        this.V_inMin = str2double(answer{2});
        this.V_inMax = str2double(answer{3});
        this.V_o = str2double(answer{4});
        this.P_o = str2double(answer{5});
        this.f_0 = str2double(answer{6});
        this.f_s = str2double(answer{7});
        this.D = str2double(answer{8});

        Converter = orderfields(this); % one-way (set) alias for cleanup
        clear this answer
        
    case 'Existing' % loads existing Converter.mat file, overwrites
        [FileName, PathName] = uigetfile('.mat');
        load([PathName, FileName])
end

clear button

%%
% *Initial Transformer Waveforms*
% 
% Waveforms from ideal transformer simulation are entered here as a starting
% point.  This is not a circuit simulator, and as such, we will leave the
% circuit sim details to a dedicated application.  Waveform data should be in
% CSV format, with columns specified in distributeWindingWaveforms.m, and each
% set should have the same number of points and sampling rate.  Further, it
% should represent a single converter switching cycle, representative of the
% converter's operation, if possible.  If there is a low-frequency component,
% one full cycle of that component should be provided, with the switching
% waveform modulating it. All primary winding waveforms should be in one file,
% while all secondary windings should be in another.
% 
% * $t$:  switching time vector in [s]
% * $dt$:  differential element of switching time in [s]
% * $t_0$:  low-frequency input time vector in [s]
% * $dt_0$:  differential element of low-frequency input time in [s]
% * $v_p(t)$:  primary winding voltage in [V]
% * $i_p(t)$:  primary winding current in [A]
% * $v_{s1}(t)$:  secondary winding 1 voltage in [V]
% * $i_{s1}(t)$:  secondary winding 1 current in [A]
% * $v_{s2}(t)$:  secondary winding 2 voltage in [V]
% * $i_{s2}(t)$:  secondary winding 2 current in [A]

this = Transformer.winding;

getString = 'Select primary winding voltage and current waveform file';
[FileName, PathName] = uigetfile('.csv', getString);
ptemp = csvread([PathName, FileName], 1, 0);
    
tVec = ptemp(:, 1) - ptemp(1, 1);

getString = 'Select secondary winding voltage and current waveform file';
[FileName, PathName] = uigetfile('.csv', getString);
stemp = csvread([PathName, FileName], 1, 1); % exclude time vector

Transformer.properties.N_wp = (size(ptemp, 2) - 1)/2;
Transformer.properties.N_ws = size(stemp, 2)/2;
Transformer.properties.N_w = Transformer.properties.N_wp + Transformer.properties.N_ws;

this.primary(1:Transformer.properties.N_wp) = struct();
this.secondary(1:Transformer.properties.N_ws) = struct();

if isequal(Converter.f_0, 0)
    [Time, this] = distributeWindingWaveforms(this, ptemp, stemp, tVec, Time, Converter.f_s);
else
    [Time, this] = distributeWindingWaveformsLF(this, ptemp, stemp, tVec, Time, Converter.f_s, Converter.f_0);
end

Transformer.winding = orderfields(this);
clear FileName PathName ptemp stemp tVec this getString

%%
% *Transformer Base Properties*
%
% These values are required to initialize several of the calculations, but
% nearly all will change as more information is obtained.  A catch-all for
% transformer properties which do not apply only to core or windings, but the
% transformer as a whole.  Flux density is estimated using voltage and iterative
% value for N to place it between 0.1 and 0.5 (common ferrite range) using a
% core of average cross-sectional area.
% 
% * $B_pk$:  peak flux density in [T]
% * $K_u$:  window utilization factor
% * $\eta$:  transformer efficiency
% * $\alpha$:  transformer regulation

this = Transformer.properties;
v = Transformer.winding.primary(1).waveform.v_p;
N = 1;
Bpk = max(abs(calculateFlux(v/N, Time.t)))/300e-6;

while Bpk > 0.5
    Bpk = max(abs(calculateFlux(v/N, Time.t)))/300e-6;
    N = N + 1;
end

this.B_pk = Bpk;
this.K_u = 0.3;
this.eta = 0.99;
this.alpha = (1 - this.eta)/2; % equate core and copper loss for max eta
Transformer.properties = orderfields(this);
clear this v N Bpk

%%
% *Core Base Parameters*
%
% NOTE:  Material selection is almost exclusively frequency-based for a given
% application.  Saturation flux density and inductance per turn are malleable,
% since they depend on shape, so we can proceed immediately and without regret.
% However, to avoid difficulty in configuration optimization, we will use
% permeability and saturation flux density as secondary factors after finding
% suitable materials based on frequency.

Bpk = Transformer.properties.B_pk;
Transformer.core.material = coreMaterial(Converter.f_s, Bpk);
Transformer.core.name = 'None';

disp('Core Material:')
disp(Transformer.core.material.main)

clear coreMaterial % prevents memory access errors
clear Bpk

%% Initial Variable Calculations
% After user entry, additional values will be needed in order to begin the
% design process.  Their calculation is handled automatically in this section.

% alias get block
thisC = Converter;
thisR = Transformer.core;
thisW = Transformer.winding;
thisP = Transformer.properties;

%%
% *Winding Calculated Parameters*
% 
% In this section, the apparent power, RMS current, and the RMS value of the
% derivative of the current are calculated.  The RMS value of the current is
% calculated by extracting any non-zero mean (the DC component), and shifting
% the AC waveform to compute I_RMS = sqrt(I_DC^2 + I_ACRMS^2).  The derivative
% of the current is approximated by a finite difference and its RMS value is
% also computed.  The RMS value of the voltage is computed in the same way.  In
% the case of low-frequency oscillation, the RMS values are computed using the
% low-frequency RMS in place of the DC value.
% 
% * $P_p$:  primary VA rating in [W]
% * $P_s$:  secondary VA rating in [W]
% * $P_t$:  total winding apparent power in [W]
% * $I_{p, RMS}$:  primary root-mean-square current in [A]
% * $I_{s, RMS}$:  secondary root-mean-square current in [A]
% * $V_{p, RMS}$:  primary root-mean-square voltage in [V]
% * $V_{s, RMS}$:  secondary root-mean-square voltage in [V]

LF = thisC.f_0 > 0;
[thisP, thisW] = analyzeWindingWaveforms(thisP, thisW, LF, Time);

%%
% *Converter Calculated Parameters*
%
% * $n$:  transformer turns ratio (matrix if multi-winding)
% * $I_o$:  nominal CCA converter output current in [A]
% * $T_s$:  converter switching period in [s]

for p = 1:thisP.N_wp
    for s = 1:thisP.N_ws
        thisC.n(p, s) = thisW.secondary(s).V_sRMS/thisW.primary(p).V_pRMS;
    end
end

thisC.I_o = thisC.P_o/thisC.V_o;
thisC.T_s = 1/thisC.f_s;

% alias set block
Converter = orderfields(thisC);
Transformer.properties = orderfields(thisP);
Transformer.core = orderfields(thisR);
Transformer.winding = orderfields(thisW);
clear thisC thisR thisW thisP

%% Design Process
% The main transformer design process begins here.

%%
% *Required Geometry Coefficient*
% 
% After this step, the core should be selected using this value and the
% Ferroxcube-provided Soft Ferrite Design Tool (SFDT).
% 
% * $K_f$:  voltage waveform factor in [cycles^(-1)]
% * $K_e$:  so-called "electrical conditions" in [W/m^5]
% * $K_g$:  core geometry coefficient in [m^5]

thisC = Converter;
thisP = Transformer.properties;
thisW = Transformer.winding;

[thisP, thisW] = computeCoreGeometry(thisC, thisP, thisW, Time);

Transformer.properties = orderfields(thisP);
Transformer.winding = orderfields(thisW);
clear thisC thisP thisW

%%
% *Run SFDT*
% 
% Cores can be checked against those available at http://www.ferroxcube.com
% on the Products page, using the Search Parameter form.  Use the following
% values in the "Transformer core selection" form, selecting the main material
% (others can be checked now or later) from the Ferrite list, and with all cores
% selected from "Core families":
% 
% * Minimum throughput power:  P_t
% * Maximum throughput power:  2*P_t (safety first)
% * Converter operating frequency:  f_s
% * Converter type:  *Half-Bridge* ~~ Push-Pull (read manual if uncertain)
% * Copper fill factor:  0.3 (0.5 is too optimistic)
% * Effective duty factor:  use D on first-pass, formula in manual thereafter
% * Ambient temperature:  25 degrees C
% * Allowed temperature rise:  min(T_C - 25, 115) (limitation of SFDT)
% * All other values at default
% 
% Click "Go" to receive a table of cores suited to the application.  Individual
% cores can be selected one at a time such that clicking "Graph Pthr(f)" will
% display a lovely plot of the throughput power and peak magnetic field
% intensity as functions of time.

thisC = Converter;
thisT = Transformer;
temp = min([thisT.core.material.main.T_c - 25, 115]);
Kg = thisT.properties.K_g*1e10;  % from m^5 to cm^5 for display
fprintf('\nSFDT Input Values:\n')
fprintf('Minimum throughput power:  %g W\n', thisT.properties.P_t)
fprintf('Maximum throughput power:  %g W\n', 2*thisT.properties.P_t)
fprintf('Converter operating frequency:  %g kHz\n', thisC.f_s*1e-3)
fprintf('Ambient temperature:  25 \260C\n')
fprintf('Allowed temperature rise:  %g \260C\n', temp)
fprintf('Converter type:  Half-Bridge --> Push-Pull\n')
fprintf('Copper fill factor:  0.3\n')
fprintf('Effective duty factor:  %g\n\n', thisC.D)
fprintf('Initial Core Geometry Coefficient :  %g cm^5\n\n', Kg)

% multipass token and D_eff calculation not yet implemented (example follows)
% if pass == 1
%     disp(fprintf('Effective duty factor:  %g', thisC.D));
% else
%     disp(fprintf('Effective duty factor:  %g', thisC.D_eff));
% end
clear temp thisC thisT Kg

%%
% *Core Selection*
% 
% This section accepts the data for the freshly-selected core, and completes
% either the initial transformer core model, or the updated core model.
% Additional data will only be needed if the initial core selection changes.
% 
% * $l_e$: effective core path length in [m]
% * $A_e$: effective core cross-sectional area in [m^2]
% * $h_c$: core window height in [m]
% * $b_c$: core window breadth in [m]
% * $h_b$: bobbin window height in [m]
% * $b_b$: bobbin window breadth in [m]
%
% From these values we can immediately calculate:
% 
% * $\mathcal{R}$: core reluctance in [H^{-1}]
% * $N_{pn}$: primary n number of turns
% * $N_{sn}$: secondary n number of turns
% * $\phi(t)$: instantaneous magnetic flux waveform
% * $B(t)$: instantaneous magnetic flux density waveform
% * $B_{pk}$: peak magnetic flux in [T]
% * $i_{m, pn}$: primary winding n magnetizing current waveform
% * $i_{m, sn}$: secondary winding n magnetizing current waveform
%
% The following values will only be available from the core's data sheet if
% paired accessories (specifically, bobbins/formers) exist for the core.  If not
% provided, these values have been approximated in selectCore.m.
% 
% * $W_a$: winding/window area in [m^2]
% * $A_p$: area product (available from data sheet if former/bobbin exists) in
% [m^4]
% * $MLT$: mean-length turn (including former/bobbin) in [m]

thisC = Transformer.core;
thisP = Transformer.properties;
thisW = Transformer.winding;

thisC = selectCore(thisC, thisP);
thisP.A_p = thisC.W_a*thisC.A_e;

fprintf('Core selected: %s\n\n', thisC.name)

nwp = thisP.N_wp;
nws = thisP.N_ws;

inputStruct = struct();
inputStruct.prompt = {'Primary winding', 'N:'; 'Secondary winding', 'N:'};
inputStruct.dlgtitle = 'Number of Turns';
inputStruct.definput = {};
inputStruct.np = nwp;
inputStruct.ns = nws;

for p = 1:nwp
    thisW.primary(p).N = ceil(thisW.primary(p).V_pRMS/(thisP.K_f*Converter.f_s*thisC.material.main.B_sat*thisC.A_e));
    inputStruct.definput{1, p} = thisW.primary(p).N;
end

for s = 1:nws
    thisW.secondary(s).N = ceil(thisW.primary(1).N*Converter.n(1, s));
    inputStruct.definput{2, s} = thisW.secondary(s).N;
end

resp = createWindingInput(inputStruct);

for p = 1:nwp
    thisW.primary(p).N = str2double(resp{p});
end

for s = 1:nws
    thisW.secondary(s).N = str2double(resp{p + s});
end

% calculate flux from primary winding voltage
% assume valid transformer config (i.e. multiple primaries have same V/turn)
if nwp > 1
    thisW.phi = zeros(size(thisW.primary(1).waveform.v_p));
    N = thisW.primary(1).N;
    thisW.phi = calculateFlux(thisW.primary(1).waveform.v_p/N, Time.t);
else
    N = thisW.primary.N;
    thisW.phi = calculateFlux(thisW.primary.waveform.v_p/N, Time.t);
end

% approximate B field in core area before windings are determined
thisW.B = thisW.phi./thisC.A_e;
thisP.B_pk = max(abs(thisW.B));
thisC.mu_r = evaluatePermeability(thisC.material.main.mu_a, thisP.B_pk);
thisC.R = thisC.l_e/(thisC.mu_r*MU_0*thisC.A_e);

% approximate primary winding magnetizing current(s)
if nwp > 1
    for idx = 1:nwp
        N = thisW.primary(idx).N;
        thisW.primary(idx).waveform.i_mp = thisW.phi.*thisC.R/N;
    end
else
    N = thisW.primary.N;
    thisW.primary.waveform.i_mp = thisW.phi.*thisC.R/N;
end

% approximate secondary winding magnetizing current(s)
if nws > 1
    for idx = 1:nws
        N = thisW.secondary(idx).N;
        thisW.secondary(idx).waveform.i_ms = thisW.phi.*thisC.R/N;
    end
else
    N = thisW.secondary.N;
    thisW.secondary.waveform.i_ms = thisW.phi.*thisC.R/N;
end

if ~thisC.W_a
    h = thisC.window.height;
    b = thisC.window.breadth;
    thisC.W_a = h*b;
end

if ~thisP.A_p
    thisP.A_p = thisC.W_a*thisC.A_e;
end

if ~thisC.MLT
    h = thisC.window.height;
    r = thisC.d_center/2; % temp radius value
    thisC.MLT = (2*pi)*(r + h/2); % mid-window circumference
end

thisP.K_g = (thisC.W_a*thisC.A_e^2*thisP.K_u)/thisC.MLT;

% display for use in LitzOpt (in millimeters)
% fprintf('LitzOpt Parameters:\n')
% fprintf('Breadth of the core window (mm):  %g\n', thisC.window.breadth*1e3)
% fprintf('Height of the core window (mm):  %g\n', thisC.window.height*1e3)
% fprintf('Breadth of the bobbin window (mm):  %g\n', thisC.bobbin.breadth*1e3)
% fprintf('Height of the bobbin window (mm):  %g\n', thisC.bobbin.height*1e3)
% fprintf('Gap length (mm):  %g\n', 0)
% fprintf('Center Leg Diameter (mm):  %g\n\n', thisC.d_center*1e3)

Transformer.core = orderfields(thisC);
Transformer.properties = orderfields(thisP);
Transformer.winding = orderfields(thisW);
clear thisC thisP thisW r h b N nwp nws selectCore p s
clear resp inputStruct

%%
% *Initial Wire Selection*
% 
% Selecting the initial wire diameter requires calculating the effective
% frequency of the waveform applied to the wire; for nonsinusoidal currents, the
% best practice is to use the ratio of the RMS value of the first derivative of
% the winding current to the RMS value of that current, all divided by 2*pi.
% This allows the calculation of the skin depth to which the current will
% permeate the wire.  Thus, to homogenize the current distribution as much as
% possible, we desire a bundle radius of no greater than 1.5 skin depths (in 
% Type I litz), and a strand diameter of about 1/4 of a skin depth.
% 
% * $f_{eff, p}$:  primary winding effective frequency in [Hz]
% * $\delta_{p}$:  primary winding skin depth in [m]
% * $d_{strand, p, max}$:  maximum primary winding strand diameter in [m]
% * $A_{strand, p, max}$:  maximum primary winding strand area in [m^2]
% * $f_{eff, s}$:  secondary winding effective frequency in [Hz]
% * $\delta_{s}$:  secondary winding skin depth in [m]
% * $d_{strand, s, max}$:  maximum secondary winding strand diameter in [m]
% * $A_{strand, s, max}$:  maximum secondary winding strand area in [m^2]
% * $t_{ins}$:  thickness of additional inter-winding space/insulation [m]

thisP = Transformer.winding.primary;
thisS = Transformer.winding.secondary;
nwp = Transformer.properties.N_wp;
nws = Transformer.properties.N_ws;

if nwp > 1
    for idx = 1:nwp
        thisP(idx).f_pEff = (1/(2*pi))*thisP(idx).di_pdt_RMS/thisP(idx).I_pRMS;
        thisP(idx).delta_p = sqrt(RHO_CU/(pi*thisP(idx).f_pEff*MU_CU));
        thisP(idx).d_pMax = thisP(idx).delta_p/4;
        thisP(idx).A_pMax = pi*thisP(idx).d_pMax^2/4;
    end
else
    thisP.f_pEff = (1/(2*pi))*thisP.di_pdt_RMS/thisP.I_pRMS;
    thisP.delta_p = sqrt(RHO_CU/(pi*thisP.f_pEff*MU_CU));
    thisP.d_pMax = thisP.delta_p/4;
    thisP.A_pMax = pi*thisP.d_pMax^2/4;    
end

if nws > 1
    for idx = 1:nws
        thisS(idx).f_sEff = (1/(2*pi))*thisS(idx).di_sdt_RMS/thisS(idx).I_sRMS;
        thisS(idx).delta_s = sqrt(RHO_CU/(pi*thisS(idx).f_sEff*MU_CU));
        thisS(idx).d_sMax = thisS(idx).delta_s/4;
        thisS(idx).A_sMax = pi*thisS(idx).d_sMax^2/4; 
    end
else
    thisS.f_sEff = (1/(2*pi))*thisS.di_sdt_RMS/thisS.I_sRMS;
    thisS.delta_s = sqrt(RHO_CU/(pi*thisS.f_sEff*MU_CU));
    thisS.d_sMax = thisS.delta_s/4;
    thisS.A_sMax = pi*thisS.d_sMax^2/4; 
end

% compute wire size guidelines and report
%TODO: Should save options...
b = Transformer.core.bobbin.breadth;
h = Transformer.core.bobbin.height;
GLstruct = guidelines(thisP, thisS, b, h);

Transformer.winding.t_ins = GLstruct.Wgap;

disp('Guidelines for Litz Wire Selection:')
Tp = table;
Ts = table;

if nwp > 1
    for idx = 1:nwp
        if idx > 1
            fprintf('\nPrimary Winding %d:\n', idx)
        else
            fprintf('Primary Winding %d:\n', idx)
        end
        
        fprintf('Current density minimum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.P(idx).AWGmin), GLstruct.P(idx).Amin*1e6)
        fprintf('Bobbin window maximum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.P(idx).AWGmax), GLstruct.P(idx).Amax*1e6)
        
        if ~isempty(GLstruct.P(idx).constructions) && ~isempty(fieldnames(GLstruct.P(idx).constructions))
            Tp = makeConsTable(Tp, GLstruct.P(idx).constructions);
        else
            fprintf('Construct custom wire arrangement for primary winding %d.\n', idx)
        end
    end
else
    fprintf('Primary Winding:\n')
    fprintf('Current density minimum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.P.AWGmin), GLstruct.P.Amin*1e6)
    fprintf('Bobbin window maximum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.P.AWGmax), GLstruct.P.Amax*1e6)
        
    if ~isempty(fieldnames(GLstruct.P.constructions))
        Tp = makeConsTable(Tp, GLstruct.P.constructions);
    else
        fprintf('Construct custom wire arrangement for primary winding.\n')
    end
end

disp(Tp)

if nws > 1
    for idx = 1:nws
        fprintf('\nSecondary Winding %d:\n', idx)
        fprintf('Current density minimum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.S(idx).AWGmin), GLstruct.S(idx).Amin*1e6)
        fprintf('Bobbin window maximum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.S(idx).AWGmax), GLstruct.S(idx).Amax*1e6)
        
        if ~isempty(GLstruct.S(idx).constructions) && ~isempty(fieldnames(GLstruct.S(idx).constructions))
            Ts = makeConsTable(Ts, GLstruct.S(idx).constructions);
        else
            fprintf('Construct custom wire arrangement for secondary winding %d.\n', idx)
        end
    end
else
    fprintf('\nSecondary Winding:\n')
    fprintf('Current density minimum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.S.AWGmin), GLstruct.S.Amin*1e6)
    fprintf('Bobbin window maximum bundle size AWG:  %d (%.2f mm^2)\n', round(GLstruct.S.AWGmax), GLstruct.S.Amax*1e6)
        
    if ~isempty(fieldnames(GLstruct.S.constructions))
        Ts = makeConsTable(Ts, GLstruct.S.constructions);
    else
        fprintf('Construct custom wire arrangement for secondary winding.\n')
    end
end

disp(Ts)

uiwait(msgbox('Click OK to continue to LitzOpt.', 'LitzOpt', 'modal'))

Transformer.winding.primary = orderfields(thisP);
Transformer.winding.secondary = orderfields(thisS);

clear thisP thisS GLstruct nwp nws b h Tp Ts

%%
% *Preparation for LitzOpt*
% In this section the time vector is reduced greatly in length to capture only
% the large changes (for which we use 'niner.m', a script that checks the second
% derivative of the supplied waveform).  The waveform with the largest number of
% changes will determine the size of all.  This is the result of setting up for
% manual entry, but now exists only because "it's always been that way".
%
% Note: in the case of low-frequency oscillation, its contribution will be very
% small in terms of high-frequency effects in the windings, so only the
% switching waveform will be used.

thisP = Transformer.winding.primary;
thisS = Transformer.winding.secondary;
nwp = Transformer.properties.N_wp;
nws = Transformer.properties.N_ws;

vectorLength = 0;
tkeep = [];
idxkeep = tkeep;

% extract and compare length of time vector for primary winding(s)
if nwp > 1
    for w = 1:nwp
        if Converter.f_0 > 0
            iwf = thisP(w).waveform.i_pHF;
            diwf = thisP(w).waveform.di_pHFdt;
            t = Time.t_s;
        else
            iwf = thisP(w).waveform.i_p;
            diwf = thisP(w).waveform.di_pdt;
            t = Time.t;
        end
        
        [time, idx] = niner(t, iwf, diwf);
        if length(time) > vectorLength
            vectorLength = length(time);
            tkeep = time;
            idxkeep = idx;
        end
    end
else
    if Converter.f_0 > 0
        iwf = thisP.waveform.i_pHF;
        diwf = thisP.waveform.di_pHFdt;
        t = Time.t_s;
    else
        iwf = thisP.waveform.i_p;
        diwf = thisP.waveform.di_pdt;
        t = Time.t;
    end
    
    [time, idx] = niner(t, iwf, diwf);
    tkeep = time;
    vectorLength = length(time);
    idxkeep = idx;
end

% extract and compare length of time vector for secondary winding(s)
if nws > 1
    for w = 1:nws
        if Converter.f_0 > 0
            iwf = thisS(w).waveform.i_sHF;
            diwf = thisS(w).waveform.di_sHFdt;
            t = Time.t_s;
        else
            iwf = thisS(w).waveform.i_s;
            diwf = thisS(w).waveform.di_sdt;
            t = Time.t;
        end
        
        [time, idx] = niner(t, iwf, diwf);
        if length(time) > vectorLength
            vectorLength = length(time);
            tkeep = time;
            idxkeep = idx;
        end
    end    
else
    if Converter.f_0 > 0
        iwf = thisS.waveform.i_sHF;
        diwf = thisS.waveform.di_sHFdt;
        t = Time.t_s;
    else
        iwf = thisS.waveform.i_s;
        diwf = thisS.waveform.di_sdt;
        t = Time.t;
    end
        
    [time, idx] = niner(t, iwf, diwf);
    if length(time) > vectorLength
        tkeep = time;
        idxkeep = idx;
    end
end

Time.t9 = tkeep;
% I = zeros(nwp + nws, length(tkeep));

if nwp > 1
    for w = 1:nwp
        if Converter.f_0 > 0
            iwf = thisP(w).waveform.i_pHF;
        else
            iwf = thisP(w).waveform.i_p;
        end
        
        thisP(w).waveform.i_p9 = iwf(idxkeep);
%         I(w, :) = thisP(w).waveform.i_p9;
    end
else
    if Converter.f_0 > 0
        iwf = thisP.waveform.i_pHF;
    else
        iwf = thisP.waveform.i_p;
    end
    
    thisP.waveform.i_p9 = iwf(idxkeep);
%     I(1, :) = thisP.waveform.i_p9;
end

if nws > 1
    for w = nwp + 1:nwp + nws
        if Converter.f_0 > 0
            iwf = thisS(w - nwp).waveform.i_sHF;
        else
            iwf = thisS(w - nwp).waveform.i_s;
        end
        
        thisS(w - nwp).waveform.i_s9 = iwf(idxkeep);
%         I(w, :) = thisS(w - nwp).waveform.i_s9;
    end
else
    if Converter.f_0 > 0
        iwf = thisS.waveform.i_sHF;
    else
        iwf = thisS.waveform.i_s;
    end
    
    thisS.waveform.i_s9 = iwf(idxkeep);
%     I(nwp + 1, :) = thisS.waveform.i_s9;
end

% tableNames = {'Duration', 'Current'};
% 
% currentTempTable = table([0, diff(Time.t9*1e6)]', I');
% currentTempTable.Properties.VariableNames = tableNames;
% 
% fprintf('\nLitzOpt Time and Current Vectors:\n');
% disp(currentTempTable)

Transformer.winding.primary = orderfields(thisP);
Transformer.winding.secondary = orderfields(thisS);
clear thisP thisS
clear time tkeep idxkeep t w idx I iwf diwf vectorLength
clear tableNames currentTempTable

%%
% *LitzOpt Optimization*
% 
% Save the file, then choose 'Existing Project' and load it again.  This calls
% 'litzoptEdit.m' and returns a table and several plots from which we can obtain
% the following values:
% 
% * $AWG_{s, pn}$:  nth primary winding strand AWG
% * $AWG_{s, sn}$:  nth secondary winding strand AWG
% * $N_{s, pn}$:  nth primary winding number of strands
% * $N_{s, sn}$: nth secondary winding number of strands
% 
% These in turn allow the calculation of:
% 
% * $d_{pn}$:  nth primary winding strand diameter in [m]
% * $d_{sn}$:  nth secondary winding 1 strand diameter in [m]
% * $A_{pn}$:  nth primary winding strand area in [m^2]
% * $A_{sn}$:  nth secondary winding strand area in [m^2]
% * $A_{Cu, pn}$:  nth primary winding copper area in [m^2]
% * $A_{Cu, sn}$:  nth secondary winding copper area in [m^2]
% * $A_{Cu}$:  total winding copper area in [m^2]
% * $AWG_{e, pn}$: nth primary winding equivalent AWG
% * $AWG_{e, sn}$: nth secondary winding equivalent AWG
%
% Each wire has a nominal outer diameter, which is approximated using a packing
% factor tuned to a known commercial product (NE wire) of 0.68125.  The primary
% winding diameter is calculated using the core dimensions and this nominal
% outer diameter.  The secondary windings are dependent on number of turns and
% outer diameter of both windings for their position information.
%TODO: compute eta in litzOptWrapper
%TODO: allow user to select one of the options

Transformer = litzOptWrapper(Transformer, Time);

% for now, enter results in the provided input dialog:
Transformer.winding = constructWinding(Transformer);

%%
% *Transformer Inductance Calculations*
% Now that the transformer is completely designed, we can evaluate its
% magnetizing inductance, mutual inductance between windings, and the leakage
% inductance in each winding.

thisW = Transformer.winding;
thisP = Transformer.properties;

% deal windings into cell array, all primary then all secondary
tempWindings = cell(1, thisP.N_w);
nwp = thisP.N_wp;
nws = thisP.N_ws;
count = 1;

if nwp > 1
    for idx = 1:nwp
        tempWindings{count} = thisW.primary(idx);
        count = count + 1;
    end
else
    tempWindings{count} = thisW.primary;
    count = count + 1;
end

if nws > 1
    for idx = 1:nws
        tempWindings{count} = thisW.secondary(idx);
        count = count + 1;
    end
else
    tempWindings{count} = thisW.secondary;
end

[thisW.M, thisW.Ml] = calculateInductance(tempWindings, Transformer.core, thisP.N_w, Converter.f_s, tol);

% assign magnetizing and leakage inductances
thisP.L_m = thisW.M(1, 1) - thisW.Ml(1, 1);

fprintf('\nInductance Calculations:\n')
fprintf('Magnetizing Inductance (H):  %g\n', thisP.L_m)

count = 1;

if nwp > 1
    for idx = 1:nwp
        thisW.primary(idx).Ll = thisW.Ml(count, count);
        count = count + 1;
        
        fprintf('Primary %d Leakage Inductance (H):  %g\n', idx, thisW.primary(idx).Ll)
    end
else
    thisW.primary.Ll = thisW.Ml(count, count);
    count = count + 1;
    
    fprintf('Primary Leakage Inductance (H):  %g\n', thisW.primary.Ll)
end

if nws > 1
    for idx = 1:nws
        thisW.secondary(idx).Ll = thisW.Ml(count, count);
        count = count + 1;
        
        fprintf('Secondary %d Leakage Inductance (H):  %g\n', idx, thisW.secondary(idx).Ll)
    end
else
    thisW.secondary.Ll = thisW.Ml(count, count);
    
    fprintf('Secondary Leakage Inductance (H):  %g\n', thisW.secondary.Ll)
end

Transformer.winding = orderfields(thisW);
Transformer.properties = orderfields(thisP);

clear thisW thisP tempWindings count nwp nws I N

%%
% *Refined Core Geometry Coefficient and Area Product*
% With the updated litz wire specification from litzOpt, Kg and Ap can be
% recalculated for a full view of the current transformer.  This allows for the
% refinement of many of the assumed values at the start and continued refinement
% of said values on subsequent iterations.  It also allows for the calculation
% of loss values for comparison.  Future work in this section could provide
% options for continued iteration, but for now manual iteration is necessary in
% order to proceed with the transformer design.

thisC = Transformer.core;
thisW = Transformer.winding;
thisP = Transformer.properties;

% find Wa required with the new litz (and check it)
thisP.W_aMin = thisW.A_Cu/thisP.K_u;
passWa = logical(thisC.W_a > thisP.W_aMin);

% recalc Kg via required Wa (and check it) 
thisP.K_gMin = thisP.W_aMin*thisC.A_e^2*thisP.K_u/thisC.MLT;
passKg = logical(thisP.K_g > thisP.K_gMin);

% calculate area product (and check it)
thisP.A_pMin = thisC.A_e*thisP.W_aMin;
passAp = logical(thisP.A_p > thisP.A_pMin);

fprintf('\nConsistency Checks:\n')
fprintf('Core has sufficient window area:  %d\n', passWa)
fprintf('Core meets K_g requirement:  %d\n', passKg)
fprintf('Core meets A_p requirement:  %d\n', passAp)

Transformer.core = orderfields(thisC);
Transformer.winding = orderfields(thisW);
Transformer.properties = orderfields(thisP);

clear thisC thisW thisP

%% Capacitance Computation
%TODO: this.

%% Loss Calculations
% In this final section of the design tool, the losses are calculated and
% presented for inspection.  Eddy effects are investigated, and copper loss is 
% calculated inclusive thereof.  coreLoss can calculate the loss due to the core
% very accurately using the iGSE method given Steinmetz parameters.  A database
% of Ferroxcube's common power conversion materials' parameters has been
% constructed for this purpose. 
% (2004 Nan/Sullivan)
% (2013 Hurley/Wolfle)
% (2014 Kazimierczuk)
%
%TODO: add optimization routine
%TODO: add state space model calculation
%TODO: add support for different types of windings and also calculate porosity
%      in windingResistance
%TODO: need to iterate over B computation as mu and L values change until 
%      delta < tol; also need to update reluctance each step

thisC = Converter;
thisR = Transformer.core;
thisW = Transformer.winding;
thisP = Transformer.properties;

nwp = thisP.N_wp;
nws = thisP.N_ws;

% final K_u calculation for this iteration
thisP.K_u = thisW.A_Cu/thisR.W_a;

% winding resistance and power loss
[thisW, thisP.P_Cu] = windingResistance(thisW);

% calculate i_m, magnetizing flux, B(t), and B_pk in core
% Note: if there is no magnetizing branch in the model, i_m will be zero
% Also: current is defined positive out of the secondary winding dotted terminal
thisW.i_m = zeros(size(Time.t));
thisR.MMF = zeros(size(Time.t));
thisR.phi = zeros(size(Time.t));

% Compute voltage across magnetizing inductance from both ends
vLmPri = zeros(size(Time.t));
vLmSec = zeros(size(Time.t));

if nwp > 1
    for idx = 1:nwp
        I = thisW.primary(idx).waveform.i_p;
        thisW.i_m = thisW.i_m + I;
        vLl = thisW.primary(idx).waveform.di_pdt*thisW.primary(idx).Ll;
        vr = thisW.primary(idx).R*I;
        vLmPri = vLmPri + (thisW.primary(idx).waveform.v_p - vLl - vr)/thisW.primary(idx).N;
    end
    
    vLmPri = vLmPri/nwp;
else
    thisW.i_m = thisW.i_m + thisW.primary.waveform.i_p;
    vLl = thisW.primary.waveform.di_pdt*thisW.primary.Ll;
    vr = thisW.primary.R*thisW.primary.waveform.i_p;
    vLmPri = vLmPri + (thisW.primary.waveform.v_p - vLl - vr)/thisW.primary.N;
end

if nws > 1
    for idx = 1:nws
        I = thisW.secondary(idx).waveform.i_s;
        N = thisW.secondary(idx).N;
        thisW.i_m = thisW.i_m - N*I;
        vLl = thisW.secondary(idx).waveform.di_sdt*thisW.secondary(idx).Ll;
        vr = thisW.secondary(idx).R*I;
        vLmSec = vLmSec + (thisW.secondary(idx).waveform.v_s + vLl + vr)/N;
    end
    
    vLmSec = vLmSec/nws;
else
    thisW.i_m = thisW.i_m - thisW.secondary.N*thisW.secondary.waveform.i_s;
    vLl = thisW.secondary.waveform.di_sdt*thisW.secondary.Ll;
    vr = thisW.secondary.R*thisW.secondary.waveform.i_s;
    vLmSec = vLmSec + (thisW.secondary.waveform.v_s + vLl + vr)/thisW.secondary.N;
end

% Compute flux, MMF, and flux density using average of voltage from primary and
% secondary windings to smooth out differences due to numerical error
thisR.phi = calculateFlux((vLmPri + vLmSec)/2, Time.t);
thisR.MMF = thisR.phi*thisR.R;
thisR.B = thisR.phi./thisR.A_e;
thisP.B_pk = max(abs(thisR.B));
thisR.mu_r = evaluatePermeability(thisR.material.main.mu_a, thisP.B_pk);

% check peak magnetic field strength against saturation value
passBsat = logical(thisP.B_pk < thisR.material.main.B_sat);
fprintf('Core meets B_sat requirement:  %d\n', passBsat)

% Window peak current density formulation based on McLyman, Hurley/Wolfle
if thisC.f_0 > 0
    thisP.J_pk = thisP.P_tHF/(thisP.K_fs*thisC.f_s*thisP.K_u*thisP.A_p*thisP.B_pk);
else
    thisP.J_pk = thisP.P_tHF/(thisP.K_f*thisC.f_s*thisP.K_u*thisP.A_p*thisP.B_pk);
end

% current density safety checks
% includes per-unit diameter surface-to-volume ratio increase 4/d
passJ = logical(thisP.J_pk < J_MAX*4);

fprintf('Window peak current density is within safety limit: %d\n', passJ)

if nwp > 1
    for idx = 1:nwp
        thisW.primary(idx).J = thisW.primary(idx).I_pRMS/thisW.primary(idx).A_Cu;
        passJ = logical(thisW.primary(idx).J < J_MAX*4);
        fprintf('Primary winding %d can support current density:  %d\n', idx, passJ)
    end
else
    thisW.primary.J = thisW.primary.I_pRMS/thisW.primary.A_Cu;
    passJ = logical(thisW.primary.J < J_MAX*4);
    fprintf('Primary winding can support current density:  %d\n', passJ)
end

if nws > 1
    for idx = 1:nws
        thisW.secondary(idx).J = thisW.secondary(idx).I_sRMS/thisW.secondary(idx).A_Cu;
        passJ = logical(thisW.secondary(idx).J < J_MAX*4);
        fprintf('Secondary winding %d can support current density:  %d\n', idx, passJ)
    end
else
    thisW.secondary.J = thisW.secondary.I_sRMS/thisW.secondary.A_Cu;
    passJ = logical(thisW.secondary.J < J_MAX*4);
    fprintf('Secondary winding can support current density:  %d\n', passJ)
end

fprintf('\nPeak Magnetic Flux Density (T):  %g\n', thisP.B_pk)
fprintf('Peak Current Density (A/mm^2):  %g\n', thisP.J_pk*1e-6)

% get Steinmetz parameters
%TODO: handle f_eff for multi-frequency situation
%TODO: model temperature dependence of k (i.e. k_T = k*(k2*T^2 + k1*T + k0))
SteinmetzOpts = thisR.material.main.Steinmetz;
freqs = [SteinmetzOpts.f];

if any(isequal(abs(freqs - thisC.f_s), 0))
    [~, opt] = min(abs(freqs - thisC.f_s));
    alpha = SteinmetzOpts(opt).alpha;
    beta = SteinmetzOpts(opt).beta;
    k = SteinmetzOpts(opt).k;
else
    alpha = SteinmetzOpts(numel(freqs)).alpha;
    beta = SteinmetzOpts(numel(freqs)).beta;
    k = SteinmetzOpts(numel(freqs)).k;
end

% assign appropriate coefficients based on B_pk and f by taking the max "loss"
if numel(alpha) > 1
    Bpk = thisP.B_pk;
    f = thisC.f_s;
    [~, idx] = max(k.*f.^alpha.*Bpk.^beta);
    alpha = alpha(idx);
    beta = beta(idx);
    k = k(idx);
end

% calculate core losses (coreloss.m will do this with iGSE given a B waveform)
% there is an error in the file that produces mW/m^3 instead of a reasonable
% value; adjustment has been made in the form of a constant factor 1e-3.
% Further amendments have been made to update the script to modern best
% practices in programming, specifically in MATLAB. Vectorization has replaced
% loops and counters, all variables have been named for self-documentation, and
% the speed has improved immensely.  The PWL approximation to k_i has also been
% replaced by the full numerical integration explained in the paper.  This
% version is included as 'corelossEdit.m'.
%
% call coreloss.m and suppress console output in favor of our formatting
thisR.P_V = corelossEdit(Time.t, thisR.B, alpha, beta, k, 1)*1e-3; % [W/m^3]
thisP.P_Fe = thisR.P_V*thisR.V_e;
thisP.P = thisP.P_Fe + thisP.P_Cu;
thisP.eta = (thisC.P_o - thisP.P)/thisC.P_o;

[thisP.T, thisP.DeltaT] = computeTempRise(thisP, thisR, thisW);

figure
plot(Time.t, thisW.B)
grid on
title('Approx. Winding Flux Density and Magnetizing Flux Density')
hold on
plot(Time.t, thisR.B)
hold off
xlabel('t [s]')
ylabel('B [T]')
legend('B_{approx.}', 'B_{core}')

% compute equivalent core resistance and display it and winding resistances
thisR.R_s = computePERMS((vLmPri + vLmSec)/2)^2/thisP.P_Fe;
fprintf('\nEquivalent Core Resistance (Ohm):  %g\n', thisR.R_s)

if nwp > 1
    for idx = 1:nwp
        fprintf('Primary %d Winding Resistance (Ohm):  %g\n', idx, thisW.primary(idx).R)
    end
else
    fprintf('Primary Winding Resistance (Ohm):  %g\n', thisW.primary.R)
end

if nws > 1
    for idx = 1:nws
        fprintf('Secondary %d Winding Resistance (Ohm):  %g\n', idx, thisW.secondary(idx).R)
    end
else
    fprintf('Secondary Winding Resistance (Ohm):  %g\n', thisW.secondary.R)
end

fprintf('\nLoss Analysis:\n')
fprintf('Copper Losses:  %g W\n', thisP.P_Cu)
fprintf('Core Losses:  %g W\n', thisP.P_Fe)
fprintf('Total Losses:  %g W\n', thisP.P)
fprintf('Transformer Efficiency:  %g %%\n', thisP.eta*100)
fprintf('Temperature Rise: %g \260C\n', thisP.DeltaT)

thisW.primary = orderfields(thisW.primary);
thisW.secondary = orderfields(thisW.secondary);

Converter = orderfields(thisC);
Transformer.core = orderfields(thisR);
Transformer.winding = orderfields(thisW);
Transformer.properties = orderfields(thisP);

% optional [select new core and provide values for Wa, Ae, le, MLT, k/alpha/beta]
% - do this when any of the 4 consistence values fail, or if the loss is too
% high, or if the change from iteration to iteration is unacceptable

clear thisC thisR thisP thisW
clear alpha beta k f freqs SteinmetzOpts N Bpk
clear I vLmPri vLmSec vLl vr nwp nws idx

%% Program Cleanup
% Here we can save the current transformer to file for future iteration or ease
% of access.  We can also clear any junk left over from the program run.
%TODO: select directory instead of saving in WD.
answer = questdlg('Save results?', 'Save or Discard', 'Yes', 'No', 'Yes');

if isequal(answer, 'Yes')
    save(['Transformer_', Transformer.core.name], 'Transformer');
    save(['Converter_', Converter.name], 'Converter');
    save(['Time_', Transformer.core.name], 'Time');
end

clear answer