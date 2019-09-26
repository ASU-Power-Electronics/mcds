%% computeCoreGeometry
% Computes Core Geometry Coefficient Kg, as well as its precursors the waveform
% factor Kf and electrical conditions Ke.  Handles both switching waveform and
% any low-frequency oscillation, and computes the total for the composite
% waveform as well.  Uses a combination of material from Hurley/Wolfle,
% Kazimierczuk, and a novel derivation of composite waveform values.

function [thisP, thisW] = computeCoreGeometry(thisC, thisP, thisW, Time)
    global SIGMA_CU
    
    if thisP.N_wp > 1
        if isequal(thisC.f_0, 0)
            vwf = thisW.primary(1).waveform.v_p;
            iwf = thisW.primary(1).waveform.i_p;
            t = Time.t;
            symmetric = any(iwf < 0);

            thisP.K_f = waveformFactor(vwf, t, symmetric);
            thisP.K_e = 0.5*SIGMA_CU*thisP.K_f^2*thisC.f_s^2*thisP.B_pk^2;
            thisP.K_g = thisP.P_t/(2*thisP.alpha*thisP.K_e); % returned to 100%
        else
            vwf = thisW.primary(1).waveform.v_pHF;
            iwf = thisW.primary(1).waveform.i_pHF;
            t = Time.t_s;
            symmetric = any(iwf < 0);

            thisP.K_fs = waveformFactor(vwf, t, symmetric);
            thisP.K_es = 0.5*SIGMA_CU*thisP.K_fs^2*thisC.f_s^2*(thisP.B_pk*thisP.P_tHF/thisP.P_t)^2;
            thisP.K_gs = thisP.P_tHF/(2*thisP.alpha*thisP.K_es);

            vwf = thisW.primary(1).waveform.v_pLF;
            iwf = thisW.primary(1).waveform.i_pLF;
            t = Time.t_0;
            symmetric = any(iwf < 0);

            thisP.K_f0 = waveformFactor(vwf, t, symmetric);
            thisP.K_e0 = 0.5*SIGMA_CU*thisP.K_f0^2*thisC.f_0^2*(thisP.B_pk*thisP.P_tLF/thisP.P_t)^2;
            thisP.K_g0 = thisP.P_tLF/(2*thisP.alpha*thisP.K_e0);

            f = sqrt(thisC.f_0*thisC.f_s);
            thisP.K_f = sqrt(thisP.K_f0*thisP.K_fs);
            thisP.K_e = 0.5*SIGMA_CU*thisP.K_f^2*f^2*thisP.B_pk^2;
            thisP.K_g = thisP.P_t/(2*thisP.K_e)*(thisP.P_tHF*thisP.K_g0*thisP.K_e0 + thisP.P_tLF*thisP.K_gs*thisP.K_es)/(thisP.P_tLF*thisP.P_tHF);
        end
    else
        if isequal(thisC.f_0, 0)
            vwf = thisW.primary.waveform.v_p;
            iwf = thisW.primary.waveform.i_p;
            t = Time.t;
            symmetric = any(iwf < 0);

            thisP.K_f = waveformFactor(vwf, t, symmetric);
            thisP.K_e = 0.5*SIGMA_CU*thisP.K_f^2*thisC.f_s^2*thisP.B_pk^2;
            thisP.K_g = thisP.P_t/(2*thisP.alpha*thisP.K_e);
        else
            vwf = thisW.primary.waveform.v_pHF;
            iwf = thisW.primary.waveform.i_pHF;
            t = Time.t_s;
            symmetric = any(iwf < 0);

            thisP.K_fs = waveformFactor(vwf, t, symmetric);
            thisP.K_es = 0.5*SIGMA_CU*thisP.K_fs^2*thisC.f_s^2*(thisP.B_pk*thisP.P_tHF/thisP.P_t)^2;
            thisP.K_gs = thisP.P_tHF/(2*thisP.alpha*thisP.K_es);

            vwf = thisW.primary.waveform.v_pLF;
            iwf = thisW.primary.waveform.i_pLF;
            t = Time.t_0;
            symmetric = any(iwf < 0);

            thisP.K_f0 = waveformFactor(vwf, t, symmetric);
            thisP.K_e0 = 0.5*SIGMA_CU*thisP.K_f0^2*thisC.f_0^2*(thisP.B_pk*thisP.P_tLF/thisP.P_t)^2;
            thisP.K_g0 = thisP.P_tLF/(2*thisP.alpha*thisP.K_e0);

            f = sqrt(thisC.f_0*thisC.f_s);
            thisP.K_f = sqrt(thisP.K_f0*thisP.K_fs);
            thisP.K_e = 0.5*SIGMA_CU*thisP.K_f^2*f^2*thisP.B_pk^2;        
            thisP.K_e = 0.5*SIGMA_CU*thisP.K_f0*thisC.f_0*thisP.K_fs*thisC.f_s*thisP.B_pk^2;
            thisP.K_g = thisP.P_t/(2*thisP.K_e)*(thisP.P_tHF*thisP.K_g0*thisP.K_e0 + thisP.P_tLF*thisP.K_gs*thisP.K_es)/(thisP.P_tLF*thisP.P_tHF);
        end
    end
end