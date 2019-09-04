%% Transformer Design Script
% *R. Scott Mongrain*
% - _June 2016-May 2019_
%
% This program is a tool for use in designing power electronics magnetic
% components (specifically transformers and inductors of various descriptions).
%
% In round figures, the design process goes as follows:
% 
% # Design converter with ideal transformer (PLECS preferred)
% # Provide design values and waveforms
% # Obtain SFDT input
% # Run SFDT, return with core values (or allow automatic selection)
% # Obtain LitzOpt input
% # Run LitzOpt, return with wire values
% # Obtain necessary modifications to Ap/Kg
% # Iterate over the previous seven steps until satisfied with delta
% # Simulate (circuit and FEM), build, test, enjoy

%% Development Notes
% * Eventually would like to connect this to converterDesign.m for a more
% complete solution to that problem.
% * Inductor design needs a solution
% * The core database is incomplete - needs more size and shape info
% * - Where did 3C98 go?
% * GUI design
% * *Transformer Structure Implementation*
% * - Can combine thisC, thisW, etc. into this as substructures (e.g. this.C)
% * *Other Stuff/Current Work*
% * - Format console output as tables (or messageboxes) where possible
% * - Implement iteration and pass counter token
% * - Implement D_eff for SFDT
% * - Implement manual overrides for windings (with UI)
% * - Generate Maxwell vbscripts based on core/windings
% * - Attach to PLECS via simulink blockset for faster iterations
% * - Investigate *webread* for scraping HTML
% * - https://blogs.mathworks.com/loren/2017/07/10/web-scraping-and-mining-unstructured-data-with-matlab/
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
J_MAX = 3e6;

%%
% *Structure Initializations*
% 
% The design documentation will describe the structures in detail.

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
% transformer.
%%
% *Converter Base Parameters*
% 
% * $V_{in, min}$:  minimum converter input voltage in [V]
% * $V_{in, max}$:  maximum converter input voltage in [V]
% * $V_{in}$: nominal converter input voltage in [V]
% * $V_o$:  nominal converter output voltage in [V]
% * $P_o$:  nominal converter output power in [W]
% * $V_f$:  converter diode forward voltage in [V]
% * $f_s$:  converter switching frequency in [Hz]
% * $D_{max}$:  maximum converter duty ratio (strict inequality)

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
                           'Switching Frequency [Hz]', ...
                           'Duty Ratio [#]'}, ...
                          'Converter Parameters', ...
                          1, ...
                          {'0', '0', '0', '0', '0', '0', '0'});
                      
        if isempty(answer) % cancel button handling
            error('Aborted, exiting design script.')
        end
        
        this = Converter; % one-way (get) alias for brevity
        this.name = answer{1};
        this.V_inMin = str2double(answer{2});
        this.V_inMax = str2double(answer{3});
        this.V_o = str2double(answer{4});
        this.P_o = str2double(answer{5});
        this.f_s = str2double(answer{6});
        this.D = str2double(answer{7});

        Converter = orderfields(this); % one-way (set) alias for cleanup
        clear this answer
        
    case 'Existing' % loads existing Converter.mat file, overwrites
        [FileName, PathName] = uigetfile('.mat');
        load([PathName, FileName])
end

clear button

%%
% *Winding Base Parameters*
% 
% * $N_w$:  number of windings
% * $N$:  transformer turns by winding
% 
%TODO: Add option to calculate number of turns automatically by
% N = ceil(V_RMS/(K_f*f*B_sat*Ae)); will likely need to defer until later

answer = inputdlg({'Number of primary windings', ...
                   'Number of turns in each primary winding (comma-separated)', ...
                   'Number of secondary windings', ...
                   'Number of turns in each secondary winding (comma-separated)'}, ...
                  'Transformer Windings and Turns', ...
                  1, ...
                  {'1', '1', '1', '1'});

thisW = Transformer.winding;
thisP = Transformer.properties;
pWindings = str2double(answer{1});
pTurns = str2num(answer{2}); % possibly non-scalar, using num instead of double
sWindings = str2double(answer{3});
sTurns = str2num(answer{4});

thisP.N_wp = pWindings;
thisP.N_ws = sWindings;
thisP.N_w = pWindings + sWindings;

if pWindings > 1
    for i = 1:pWindings
        thisW.primary(i).N = pTurns(i);
    end
else
    thisW.primary.N = pTurns;
end

if sWindings > 1
    for i = 1:sWindings
        thisW.secondary(i).N = sTurns(i);
    end
else
    thisW.secondary.N = sTurns;
end

Transformer.winding = orderfields(thisW);
Transformer.properties = orderfields(thisP);
clear thisW thisP answer pWindings pTurns sWindings sTurns

%%
% *Initial Transformer Waveforms*
% 
% Waveforms from ideal transformer simulation are entered here as a starting
% point.  This is not a circuit simulator, and as such, we will leave the
% circuit sim details to a dedicated application.  Waveform data should be in
% CSV format, with columns specified in distributeWindingWaveforms.m, and each
% set should have the same number of points and sampling rate.  Further, it
% should represent a single converter switching cycle, representative of the
% converter's operation, if possible.  All primary winding waveforms should be
% in one file, while all secondary windings should be in another.
% 
% * $t$:  time vector in [s]
% * $dt$:  differential element of time in [s]
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
Time.t = linspace(0, 1/Converter.f_s, 1025);
Time.dt = Time.t(2) - Time.t(1);

getString = 'Select secondary winding voltage and current waveform file';
[FileName, PathName] = uigetfile('.csv', getString);
stemp = csvread([PathName, FileName], 1, 1); % exclude time vector for brevity

this = distributeWindingWaveforms(this, ptemp, stemp, tVec, Time.t);

Transformer.winding = orderfields(this);
clear FileName PathName ptemp stemp tVec this getString

%%
% *Core Base Parameters*
%
% NOTE:  Material selection is almost exclusively frequency-based for a given
% application, so we can proceed immediately and without regret.

Transformer.core.material = coreMaterial(Converter.f_s);
Transformer.core.name = 'None';

disp('Core Material:')
disp(Transformer.core.material.main)

clear coreMaterial % prevents memory access errors

%%
% *Transformer Base Properties*
%
% These values are required to initialize several of the calculations, but
% nearly all will change as more information is obtained.  A catch-all for
% transformer properties which do not apply only to core or windings, but the
% transformer as a whole.
% 
% * $B_pk$:  peak flux density in [T]
% * $K_u$:  window utilization factor
% * $\eta$:  transformer efficiency
% * $\alpha$:  transformer regulation

this = Transformer.properties;
this.B_pk = 0.15;
this.K_u = 0.3;
this.eta = 0.99;
this.alpha = (1 - this.eta)/2; % equate core and copper loss for max eta
Transformer.properties = orderfields(this);
clear this

%% Initial Variable Calculations
% After user entry, additional values will be needed in order to begin the
% design process.  Their calculation is handled automatically in this section.

% alias get block
thisC = Converter;
thisR = Transformer.core;
thisW = Transformer.winding;
thisP = Transformer.properties;

%%
% *Converter Calculated Parameters*
%
% * $n$:  transformer turns ratio
% * $I_o$:  nominal converter output current in [A]
% * $V_{in}$:  nominal converter input voltage in [V]
% * $T_s$:  converter switching period in [s]
% * $D_{min}$:  minimum converter duty ratio
% * $D$:  converter duty ratio

%TODO: move to converter selection as type drop-down and add other isolated
%converter types (or remove entirely...)

% Asymmetric Half-Bridge (comments include symmetric HB)
% thisC.V_f = 1.4; % diode forward voltage drop (cut-in model)

if numel([thisW.secondary(:).N]) > 1
    thisC.n = thisW.secondary(1).N/thisW.primary.N;
else
    thisC.n = thisW.secondary.N/thisW.primary.N;
end

thisC.I_o = thisC.P_o/thisC.V_o;
% thisC.V_in = (thisC.V_inMin + thisC.V_inMax)/2;
% thisC.V_in = 34.6;
thisC.T_s = 1/thisC.f_s;
% thisC.D_min = 1/2 - sqrt(1/4 - (thisC.V_o + thisC.V_f)/(2*thisC.n*thisC.V_inMax));
% thisC.D = 1/2 - sqrt(1/4 - (thisC.V_o + thisC.V_f)/(2*thisC.n*thisC.V_in));
% thisC.D_min = (thisC.V_o + thisC.V_f)/(thisC.n*thisC.V_inMax); % symmetric
% thisC.D = thisC.V_o/(thisC.n*thisC.V_in); % symmetric
% thisC.D = 0.321;                                                                 % MANUAL EDIT

%%
% *Winding Calculated Parameters*
% 
% In this section, the apparent power, RMS current, and the RMS value of the
% derivative of the current are calculated.  The RMS value of the current is
% calculated by extracting any non-zero mean (the DC component), and shifting
% the AC waveform to compute I_RMS = sqrt(I_DC^2 + I_ACRMS^2).  The derivative
% of the current is approximated by a finite difference and its RMS value is
% computed using MATLAB's built-in rms() function.
% 
% * $P_p$:  primary VA rating in [W]
% * $P_s$:  secondary VA rating in [W]
% * $P_t$:  total winding apparent power in [W]
% * $I_{p, RMS}$:  primary root-mean-square current in [A]
% * $I_{s, RMS}$:  secondary root-mean-square current in [A]

nwp = thisP.N_wp;
nws = thisP.N_ws;
thisP.P_t = 0; % initialize throughput power to zero for summation

if nwp > 1
    for i = 1:nwp
        iwf = thisW.primary(i).waveform.i_p;
        vwf = thisW.primary(i).waveform.v_p;
        thisW.primary(i).I_pRMS = computePERMS(iwf);
        thisW.primary(i).V_pRMS = computePERMS(vwf);
        thisW.primary(i).VA = thisW.primary(i).I_pRMS*thisW.primary(i).V_pRMS;
        thisP.P_t = thisP.P_t + thisW.primary(i).VA;
        temp = [thisW.primary(i).waveform.i_p, thisW.primary(i).waveform.i_p(1)];
        thisW.primary(i).waveform.di_pdt = diff(temp)/Time.dt;
        thisW.primary(i).di_pdt_RMS = rms(thisW.primary(i).waveform.di_pdt(1:end - 1));
    end
else
    iwf = thisW.primary.waveform.i_p;
    vwf = thisW.primary.waveform.v_p;
    thisW.primary.I_pRMS = computePERMS(iwf);
    thisW.primary.V_pRMS = computePERMS(vwf);
    thisW.primary.VA = thisW.primary.I_pRMS*thisW.primary.V_pRMS;
    thisP.P_t = thisP.P_t + thisW.primary.VA;
    temp = [thisW.primary.waveform.i_p, thisW.primary.waveform.i_p(1)];
    thisW.primary.waveform.di_pdt = diff(temp)/Time.dt;
    thisW.primary.di_pdt_RMS = rms(thisW.primary.waveform.di_pdt(1:end - 1));
end

if nws > 1
    for i = 1:nws
        iwf = thisW.secondary(i).waveform.i_s;
        vwf = thisW.secondary(i).waveform.v_s;
        thisW.secondary(i).I_sRMS = computePERMS(iwf);
        thisW.secondary(i).V_sRMS = computePERMS(vwf);
        thisW.secondary(i).VA = thisW.secondary(i).I_sRMS*thisW.secondary(i).V_sRMS;
        thisP.P_t = thisP.P_t + thisW.secondary(i).VA;
        temp = [thisW.secondary(i).waveform.i_s, thisW.secondary(i).waveform.i_s(1)];
        thisW.secondary(i).waveform.di_sdt = diff(temp)/Time.dt;
        thisW.secondary(i).di_sdt_RMS = rms(thisW.secondary(i).waveform.di_sdt(1:end - 1));
    end
else
    iwf = thisW.secondary.waveform.i_s;
    vwf = thisW.secondary.waveform.v_s;
    thisW.secondary.I_sRMS = computePERMS(iwf);
    thisW.secondary.V_sRMS = computePERMS(vwf);
    thisW.secondary.VA = thisW.secondary.I_sRMS*thisW.secondary.V_sRMS;
    thisP.P_t = thisP.P_t + thisW.secondary.VA;
    temp = [thisW.secondary.waveform.i_s, thisW.secondary.waveform.i_s(1)];
    thisW.secondary.waveform.di_sdt = diff(temp)/Time.dt;
    thisW.secondary.di_sdt_RMS = rms(thisW.secondary.waveform.di_sdt(1:end - 1));
end

% alias set block
Converter = orderfields(thisC);
Transformer.properties = orderfields(thisP);
Transformer.core = orderfields(thisR);
Transformer.winding = orderfields(thisW);
clear thisC thisR thisW thisP nwp nws iwf vwf temp

%% Design Process
% A separate PDF will detail the design process.

%%
% *Required Geometry Coefficient*
% 
% After this step, the core should be selected using this value and the
% Ferroxcube-provided Soft Ferrite Design Tool (SFDT).
% 
% * $K_f$:  voltage waveform factor in [cycles^(-1)]
% * $K_e$:  so-called "electrical conditions" in [W/m^5]
% * $K_g$:  core geometry coefficient in [m^5]

this = Transformer.properties;

if this.N_wp > 1
    nwp = Transformer.winding.primary(1).N;
    vwf = Transformer.winding.primary(1).waveform.v_p/nwp;
    iwf = Transformer.winding.primary(1).waveform.i_p;
else
    nwp = Transformer.winding.primary.N;
    vwf = Transformer.winding.primary.waveform.v_p/nwp;
    iwf = Transformer.winding.primary.waveform.i_p;
end

if any(iwf < 0)
    symmetric = true;
else
    symmetric = false;
end

% this.K_f = 2*sqrt(2/Converter.D); % Half-bridge specific
this.K_f = waveformFactor(vwf, Time.t, symmetric);
this.K_e = 0.5*SIGMA_CU*this.K_f^2*Converter.f_s^2*this.B_pk^2;
this.K_g = this.P_t/(2*this.alpha*this.K_e)*100; % returned to 100%
Transformer.properties = orderfields(this);
clear this vwf iwf nwp symmetric

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
% * Maximum throughput power:  sqrt(2)*P_t <-- *Half-bridge*
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
Kg = thisT.properties.K_g*1e15;  % from m^5 to mm^5 for display
fprintf('\nSFDT Input Values:\n')
fprintf('Minimum throughput power:  %g W\n', thisT.properties.P_t)
fprintf('Maximum throughput power:  %g W\n', sqrt(2)*thisT.properties.P_t)
fprintf('Converter operating frequency:  %g kHz\n', thisC.f_s*1e-3)
fprintf('Ambient temperature:  25 degrees C\n')
fprintf('Allowed temperature rise:  %g degrees C\n', temp)
fprintf('Converter type:  Half-Bridge --> Push-Pull\n')
fprintf('Copper fill factor:  0.3\n')
fprintf('Effective duty factor:  %g\n\n', thisC.D)
fprintf('Initial Core Geometry Coefficient :  %g mm^5\n\n', Kg)

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
thisC.R = thisC.l_e/(thisC.material.main.mu_i*MU_0*thisC.A_e);
thisP.A_p = thisC.W_a*thisC.A_e;

fprintf('Core selected: %s\n\n', thisC.name)

nwp = thisP.N_wp;
nws = thisP.N_ws;

% calculate flux from primary winding voltage(s)
if nwp > 1
    thisW.phi = zeros(size(thisW.primary(1).waveform.v_p));
    for idx = 1:nwp
        N = thisW.primary(idx).N;
        thisW.phi = thisW.phi + calculateFlux(thisW.primary(idx).waveform.v_p/N, Time.t);
    end
else
    N = thisW.primary.N;
    thisW.phi = calculateFlux(thisW.primary.waveform.v_p/N, Time.t);
end

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

% approximate B field in core area before windings are determined
thisW.B = thisW.phi./thisC.A_e;
thisP.B_pk = max(abs(thisW.B));

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
fprintf('LitzOpt Parameters:\n')
fprintf('Breadth of the core window (mm):  %g\n', thisC.window.breadth*1e3)
fprintf('Height of the core window (mm):  %g\n', thisC.window.height*1e3)
fprintf('Breadth of the bobbin window (mm):  %g\n', thisC.bobbin.breadth*1e3)
fprintf('Height of the bobbin window (mm):  %g\n', thisC.bobbin.height*1e3)
fprintf('Gap length (mm):  %g\n', 0)
fprintf('Center Leg Diameter (mm):  %g\n\n', thisC.d_center*1e3)

Transformer.core = orderfields(thisC);
Transformer.properties = orderfields(thisP);
Transformer.winding = orderfields(thisW);
clear thisC thisP thisW r h b N nwp nws selectCore

%%
% *Initial Wire Selection*
% 
% Selecting the initial wire diameter requires calculating the effective
% frequency of the waveform applied to the wire; for nonsinusoidal currents, the
% best practice is to use the ratio of the RMS value of the first derivative of
% the winding current to the RMS value of that current, all divided by 2*pi.
% This allows the calculation of the skin depth to which the current will
% permeate the wire.  Thus, to homogenize the current distribution as much as
% possible, we desire a conductor radius of no greater than 1.5 skin depths.
% 
% * $f_{eff, p}$:  primary winding effective frequency in [Hz]
% * $\delta_{p}$:  primary winding skin depth in [m]
% * $d_{strand, p, max}$:  maximum primary winding strand diameter in [m]
% * $A_{strand, p, max}$:  maximum primary winding strand area in [m^2]
% * $f_{eff, s}$:  secondary winding effective frequency in [Hz]
% * $\delta_{s}$:  secondary winding skin depth in [m]
% * $d_{strand, s, max}$:  maximum secondary winding strand diameter in [m]
% * $A_{strand, s, max}$:  maximum secondary winding strand area in [m^2]

thisP = Transformer.winding.primary;
thisS = Transformer.winding.secondary;
nwp = Transformer.properties.N_wp;
nws = Transformer.properties.N_ws;

if nwp > 1
    for idx = 1:nwp
        thisP(idx).f_pEff = (1/(2*pi))*thisP(idx).di_pdt_RMS/thisP(idx).I_pRMS;
        thisP(idx).delta_p = sqrt(RHO_CU/(pi*thisP(idx).f_pEff*MU_CU));
        thisP(idx).d_pMax = 3*thisP(idx).delta_p;
        thisP(idx).A_pMax = pi*thisP(idx).d_pMax^2/4;
    end
else
    thisP.f_pEff = (1/(2*pi))*thisP.di_pdt_RMS/thisP.I_pRMS;
    thisP.delta_p = sqrt(RHO_CU/(pi*thisP.f_pEff*MU_CU));
    thisP.d_pMax = 2*thisP.delta_p;
    thisP.A_pMax = pi*thisP.d_pMax^2/4;    
end

if nws > 1
    for idx = 1:nws
        thisS(idx).f_sEff = (1/(2*pi))*thisS(idx).di_sdt_RMS/thisS(idx).I_sRMS;
        thisS(idx).delta_s = sqrt(RHO_CU/(pi*thisS(idx).f_sEff*MU_CU));
        thisS(idx).d_sMax = 3*thisS(idx).delta_s;
        thisS(idx).A_sMax = pi*thisS(idx).d_sMax^2/4; 
    end
else
    thisS.f_sEff = (1/(2*pi))*thisS.di_sdt_RMS/thisS.I_sRMS;
    thisS.delta_s = sqrt(RHO_CU/(pi*thisS.f_sEff*MU_CU));
    thisS.d_sMax = 2*thisS.delta_s;
    thisS.A_sMax = pi*thisS.d_sMax^2/4; 
end

% compute wire size guidelines and report
%TODO: Fix this to actually suggest something based on current density.
b = Transformer.core.bobbin.breadth*1e3; % to [mm]
GLstruct = guidelines(thisP, thisS, b);

disp('Guidelines for Litz Wire Selection:')

if nwp > 1
    for idx = 1:nwp
        fprintf('Primary Winding %d:\n', idx)
        fprintf('Safety limit for conductor bundle equivalent AWG:  %d\n', GLstruct.P(idx).minbAWG)
        fprintf('Suggested strand size AWG:  %d\n', GLstruct.P(idx).minAWG)
        fprintf('Associated number of strands:  %d\n', GLstruct.P(idx).ndMax)
        fprintf('Resulting outer diameter:  %f mm\n', GLstruct.P(idx).Db*1e3)
        fprintf('Fits %d conductor bundles per layer.\n\n', GLstruct.P(idx).nbpl)
    end
else
    fprintf('Primary Winding:\n')
    fprintf('Safety limit for conductor bundle equivalent AWG:  %d\n', GLstruct.P.minbAWG)
    fprintf('Suggested strand size AWG:  %d\n', GLstruct.P.minAWG)
    fprintf('Associated number of strands:  %d\n', GLstruct.P.ndMax)
    fprintf('Resulting outer diameter:  %f mm\n', GLstruct.P.Db*1e3)
    fprintf('Fits %d conductor bundles per layer.\n\n', GLstruct.P.nbpl)
end

if nws > 1
    for idx = 1:nws
        fprintf('Secondary Winding %d:\n', idx)
        fprintf('Safety limit for conductor bundle equivalent AWG:  %d\n', GLstruct.S(idx).minbAWG)
        fprintf('Suggested strand size AWG:  %d\n', GLstruct.S(idx).minAWG)
        fprintf('Associated number of strands:  %d\n', GLstruct.S(idx).ndMax)
        fprintf('Resulting outer diameter:  %f mm\n', GLstruct.S(idx).Db*1e3)
        fprintf('Fits %d conductor bundles per layer.\n\n', GLstruct.S(idx).nbpl)
    end
else
    fprintf('Secondary Winding:\n')
    fprintf('Safety limit for conductor bundle equivalent AWG:  %d\n', GLstruct.S.minbAWG)
    fprintf('Suggested strand size AWG:  %d\n', GLstruct.S.minAWG)
    fprintf('Associated number of strands:  %d\n', GLstruct.S.ndMax)
    fprintf('Resulting outer diameter:  %f mm\n', GLstruct.S.Db*1e3)
    fprintf('Fits %d conductor bundles per layer.\n\n', GLstruct.S.nbpl)
end

Transformer.winding.primary = orderfields(thisP);
Transformer.winding.secondary = orderfields(thisS);

clear thisP thisS GLstruct nwp nws b

%%
% *Preparation for LitzOpt*
% In this section the time vector is reduced greatly in length to capture only
% the large changes (for which we use 'niner.m', a script that checks the second
% derivative of the supplied waveform).  The waveform with the largest number of
% changes will determine the size of all.  This greatly facilitates entry of the
% data points into LitzOpt, which requires each duration and current value to be
% entered manually.  This will later be adjusted such that the file and
% structure system common to this program are usable with LitzOpt directly.  An
% added effect will be that waveforms of more than a handful of points can be
% used as they will be automatically entered.

thisP = Transformer.winding.primary;
thisS = Transformer.winding.secondary;
nwp = Transformer.properties.N_wp;
nws = Transformer.properties.N_ws;

vectorLength = 0;
tkeep = [];
idxkeep = tkeep;

% extract and compare length of time vector for primary winding(s)
if nwp > 1
    for i = 1:nwp
        [time, idx] = niner(Time.t, thisP(i).waveform.i_p, thisP(i).waveform.di_pdt);
        if length(time) > vectorLength
            vectorLength = length(time);
            tkeep = time;
            idxkeep = idx;
        end
    end
else
    [time, idx] = niner(Time.t, thisP.waveform.i_p, thisP.waveform.di_pdt);
    tkeep = time;
    vectorLength = length(time);
    idxkeep = idx;
end

% extract and compare length of time vector for secondary winding(s)
if nws > 1
    for i = 1:nws
        [time, idx] = niner(Time.t, thisS(i).waveform.i_s, thisS(i).waveform.di_sdt);
        if length(time) > vectorLength
            vectorLength = length(time);
            tkeep = time;
            idxkeep = idx;
        end
    end    
else
    [time, idx] = niner(Time.t, thisS(1).waveform.i_s, thisS(1).waveform.di_sdt);
    if length(time) > vectorLength
        tkeep = time;
        idxkeep = idx;
    end
end

Time.t9 = tkeep;
I = zeros(nwp + nws, length(tkeep));

if nwp > 1
    for i = 1:nwp
        thisP(i).waveform.i_p9 = thisP(i).waveform.i_p(idxkeep);
        I(i, :) = thisP(i).waveform.i_p9;
    end
else
    thisP.waveform.i_p9 = thisP.waveform.i_p(idxkeep);
    I(1, :) = thisP.waveform.i_p9;
end

if nws > 1
    for i = nwp + 1:nwp + nws
        thisS(i - nwp).waveform.i_s9 = thisS(i - nwp).waveform.i_s(idxkeep);
        I(i, :) = thisS(i - nwp).waveform.i_s9;
    end
else
    thisS.waveform.i_s9 = thisS.waveform.i_s(idxkeep);
    I(nwp + 1, :) = thisS.waveform.i_s9;
end

tableNames = {'Duration', 'Current'};

currentTempTable = table([0, diff(Time.t9*1e6)]', I');
currentTempTable.Properties.VariableNames = tableNames;

disp('LitzOpt Time and Current Vectors:')
disp(currentTempTable)

Transformer.winding.primary = orderfields(thisP);
Transformer.winding.secondary = orderfields(thisS);
clear thisP thisS time tkeep idxkeep idx tableNames currentTempTable I vectorLength

%%
% *LitzOpt Optimization*
% (aka future site of litzOpt call and return)
% 
% Follow the instructions in the web-based tool for now.  On the first pass,
% select "2-D Internal Field Sim", "Assume Layered Windings", enter the number
% of windings (for the selected example, 3), and choose "Piecewise linear" for
% the "Current Waveform Specification Mode".  Enter the number of time steps in
% the text field and click "Go!".  All information required has been presented,
% and no unit conversions are necessary; the values are converted before
% display.  When the values provided have been entered, click "Go!" at the
% bottom without changing any other values.  This returns a table from which we
% can obtain the following values:
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

%TODO: implement call to LitzOpt
%TODO: in the interim read from file instead of manual entry in web version
% for now, enter results here:

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

[thisW.M, thisW.Ml] = calculateInductance(tempWindings, Transformer.core, thisP.N_w, Converter.f_s);

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
% of loss values for comparison.  Future work in this section will provide
% options for continued iteration, but for now manual iteration is necessary in
% order to proceed with the transformer design.

thisC = Transformer.core;
thisW = Transformer.winding;
thisP = Transformer.properties;

% find Wa required with the new litz (and check it)
thisP.W_aMin = thisW.A_Cu/thisP.K_u;
passWa = logical(thisC.W_a > thisP.W_aMin);

% using litzOpt sizes, recalc Kg via required Wa (and check it) 
thisP.K_gMin = thisP.W_aMin*thisC.A_e^2*thisP.K_u/thisC.MLT;
passKg = logical(thisP.K_g > thisP.K_gMin);

% calculate area product (and check it)
thisP.A_pMin = thisC.A_e*thisP.W_aMin;
passAp = logical(thisP.A_p > thisP.A_pMin);

% check peak magnetic field strength against saturation value
passBsat = logical(thisP.B_pk < thisC.material.main.B_sat);

fprintf('\nConsistency Checks:\n')
fprintf('Core has sufficient window area:  %d\n', passWa)
fprintf('Core meets K_g requirement:  %d\n', passKg)
fprintf('Core meets A_p requirement:  %d\n', passAp)
fprintf('Core meets B_sat requirement:  %d\n', passBsat)

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
% Note: if there is no magnetizing branch, i_m will be zero
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
        vLmPri = vLmPri + thisW.primary(Idx).waveform.v_p - vLl - vr;
    end
else
    thisW.i_m = thisW.i_m + thisW.primary.waveform.i_p;
    vLl = thisW.primary.waveform.di_pdt*thisW.primary.Ll;
    vr = thisW.primary.R*thisW.primary.waveform.i_p;
    vLmPri = vLmPri + thisW.primary.waveform.v_p - vLl - vr;
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

% formulation based on McLyman, Hurley/Wolfle
thisP.J_pk = thisP.P_t/(thisP.K_f*thisC.f_s*thisP.K_u*thisP.A_p*thisP.B_pk);

% current density safety checks
passJ = logical(thisP.J_pk < J_MAX); % window current density

fprintf('Window peak current density is within safety limit: %d\n', passJ)

if nwp > 1
    for idx = 1:nwp
        thisW.primary(idx).J = thisW.primary(idx).I_pRMS/thisW.primary(idx).A_Cu;
        passJ = logical(thisW.primary(idx).J < J_MAX);
        fprintf('Primary winding %d can support current density:  %d\n', idx, passJ)
    end
else
    thisW.primary.J = thisW.primary.I_pRMS/thisW.primary.A_Cu;
    passJ = logical(thisW.primary.J < J_MAX);
    fprintf('Primary winding can support current density:  %d\n', passJ)
end

if nws > 1
    for idx = 1:nws
        thisW.secondary(idx).J = thisW.secondary(idx).I_sRMS/thisW.secondary(idx).A_Cu;
        passJ = logical(thisW.secondary(idx).J < J_MAX);
        fprintf('Secondary winding %d can support current density:  %d\n', idx, passJ)
    end
else
    thisW.secondary.J = thisW.secondary.I_sRMS/thisW.secondary.A_Cu;
    passJ = logical(thisW.secondary.J < J_MAX);
    fprintf('Secondary winding can support current density:  %d\n', passJ)
end

fprintf('\nPeak Magnetic Flux Density (T):  %g\n', thisP.B_pk)
fprintf('Peak Current Density (A/mm^2):  %g\n', thisP.J_pk*1e-6)

% get Steinmetz parameters (approximate linearly if needed...)
SteinmetzOpts = thisR.material.main.Steinmetz;
freqs = [SteinmetzOpts.f];
[~, idx] = min(abs(freqs - thisC.f_s));

if numel(idx) > 1
    alpha = mean([SteinmetzOpts(idx).alpha]);
    beta = mean([SteinmetzOpts(idx).beta]);
    k = mean([SteinmetzOpts(idx).k]);
else
    alpha = SteinmetzOpts(idx).alpha;
    beta = SteinmetzOpts(idx).beta;
    k = SteinmetzOpts(idx).k;
end

% calculate core losses (coreloss.m will do this with iGSE given a B waveform)
% there is an error in the file that produces mW/m^3 instead of a reasonable
% value; adjustment has been made in the form of a constant factor 1e-3.
% call coreloss.m and suppress console output in favor of our formatting
thisR.P_V = coreloss(Time.t, thisR.B, alpha, beta, k, 1)*1e-3; % [W/m^3]
thisP.P_Fe = thisR.P_V*thisR.V_e;
thisP.P = thisP.P_Fe + thisP.P_Cu;
thisP.eta = (thisC.P_o - thisP.P)/thisC.P_o;

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
thisR.R_s = rms((vLmPri + vLmSec)/2)^2/thisP.P_Fe;
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
fprintf('Copper Losses (W):  %g\n', thisP.P_Cu)
fprintf('Core Losses (W):  %g\n', thisP.P_Fe)
fprintf('Total Losses (W):  %g\n', thisP.P)
fprintf('Transformer Efficiency:  %g %%\n', thisP.eta*100)

thisW.primary = orderfields(thisW.primary);
thisW.secondary = orderfields(thisW.secondary);

Converter = orderfields(thisC);
Transformer.core = orderfields(thisR);
Transformer.winding = orderfields(thisW);
Transformer.properties = orderfields(thisP);

% optional [select new core and provide values for Wa, Ae, le, MLT, k/alpha/beta]
% - do this when any of the 4 consistence values fail, or if the loss is too
% high, or if the change from iteration to iteration is unacceptable

clear thisC thisR thisP thisW alpha beta k freqs vLmPri vLmSec vLl vr nwp nws idx SteinmetzOpts

%% Program Cleanup
% Here we can save the current transformer to file for future iteration or ease
% of access.  We can also clear any junk left over from the program run.
answer = questdlg('Save results?', 'Save or Discard', 'Yes', 'No', 'Yes');

if isequal(answer, 'Yes')
    save(['Transformer_', Transformer.core.name], 'Transformer');
    save(['Converter_', Converter.name], 'Converter');
    save(['Time_', Transformer.core.name], 'Time');
end

clear answer