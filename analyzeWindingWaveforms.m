%% analyzeWindingWaveforms
% Performs various analyses of waveforms in windings, including determining
% appropriate RMS quantities, apparent power, and peak values.  Handles pure
% switching waveforms with or without DC component, or composite waveforms with
% low-frequency modulation/carrier.
function [P, W] = analyzeWindingWaveforms(P, W, LF, T)
    nwp = P.N_wp;
    nws = P.N_ws;
    
    % initialize throughput power(s) to zero for summation
    P.P_t = 0;
    P.P_tLF = 0;
    P.P_tHF = 0;

    if ~LF
        if nwp > 1
            for p = 1:nwp
                iwf = W.primary(p).waveform.i_p;
                vwf = W.primary(p).waveform.v_p;
                W.primary(p).I_pPk = max(abs(iwf));
                W.primary(p).I_pRMS = computePERMS(iwf);
                W.primary(p).V_pRMS = computePERMS(vwf);
                W.primary(p).VA = W.primary(p).I_pRMS*W.primary(p).V_pRMS;
                P.P_t = P.P_t + W.primary(p).VA;
                temp = [W.primary(p).waveform.i_p, W.primary(p).waveform.i_p(1)];
                W.primary(p).waveform.di_pdt = diff(temp)/T.dt;
                W.primary(p).di_pdt_RMS = computePERMS(W.primary(p).waveform.di_pdt(1:end - 1));
            end
        else
            iwf = W.primary.waveform.i_p;
            vwf = W.primary.waveform.v_p;
            W.primary.I_pPk = max(abs(iwf));
            W.primary.I_pRMS = computePERMS(iwf);
            W.primary.V_pRMS = computePERMS(vwf);
            W.primary.VA = W.primary.I_pRMS*W.primary.V_pRMS;
            P.P_t = P.P_t + W.primary.VA;
            temp = [W.primary.waveform.i_p, W.primary.waveform.i_p(1)];
            W.primary.waveform.di_pdt = diff(temp)/T.dt;
            W.primary.di_pdt_RMS = computePERMS(W.primary.waveform.di_pdt(1:end - 1));
        end

        if nws > 1
            for s = 1:nws
                iwf = W.secondary(s).waveform.i_s;
                vwf = W.secondary(s).waveform.v_s;
                W.secondary(s).I_sPk = max(abs(iwf));
                W.secondary(s).I_sRMS = computePERMS(iwf);
                W.secondary(s).V_sRMS = computePERMS(vwf);
                W.secondary(s).VA = W.secondary(s).I_sRMS*W.secondary(s).V_sRMS;
                P.P_t = P.P_t + W.secondary(s).VA;
                temp = [W.secondary(s).waveform.i_s, W.secondary(s).waveform.i_s(1)];
                W.secondary(s).waveform.di_sdt = diff(temp)/T.dt;
                W.secondary(s).di_sdt_RMS = computePERMS(W.secondary(s).waveform.di_sdt(1:end - 1));
            end
        else
            iwf = W.secondary.waveform.i_s;
            vwf = W.secondary.waveform.v_s;
            W.secondary.I_sPk = max(abs(iwf));
            W.secondary.I_sRMS = computePERMS(iwf);
            W.secondary.V_sRMS = computePERMS(vwf);
            W.secondary.VA = W.secondary.I_sRMS*W.secondary.V_sRMS;
            P.P_t = P.P_t + W.secondary.VA;
            temp = [W.secondary.waveform.i_s, W.secondary.waveform.i_s(1)];
            W.secondary.waveform.di_sdt = diff(temp)/T.dt;
            W.secondary.di_sdt_RMS = computePERMS(W.secondary.waveform.di_sdt(1:end - 1));
        end
        
        P.P_tHF = P.P_t;
    else
        if nwp > 1
            for p = 1:nwp
                % LF
                iwf = W.primary(p).waveform.i_pLF;
                vwf = W.primary(p).waveform.v_pLF;
                W.primary(p).I_pLFPk = max(abs(iwf));
                W.primary(p).I_pLFRMS = computePERMS(iwf);
                W.primary(p).V_pLFRMS = computePERMS(vwf);
                W.primary(p).VALF = W.primary(p).I_pLFRMS*W.primary(p).V_pLFRMS;
                P.P_tLF = P.P_tLF + W.primary(p).VALF;
                temp = [W.primary(p).waveform.i_pLF, W.primary(p).waveform.i_pLF(1)];
                W.primary(p).waveform.di_pLFdt = diff(temp)/T.dt_0;
                W.primary(p).di_pLFdt_RMS = computePERMS(W.primary(p).waveform.di_pLFdt(1:end - 1));
                
                % HF
                iwf = W.primary(p).waveform.i_pHF;
                vwf = W.primary(p).waveform.v_pHF;
                W.primary(p).I_pHFPk = max(abs(iwf));
                W.primary(p).I_pHFRMS = computePERMS(iwf);
                W.primary(p).V_pHFRMS = computePERMS(vwf);
                W.primary(p).VAHF = W.primary(p).I_pHFRMS*W.primary(p).V_pHFRMS;
                P.P_tHF = P.P_tHF + W.primary(p).VAHF;
                temp = [W.primary(p).waveform.i_pHF, W.primary(p).waveform.i_pHF(1)];
                W.primary(p).waveform.di_pHFdt = diff(temp)/T.dt_s;
                W.primary(p).di_pHFdt_RMS = computePERMS(W.primary(p).waveform.di_pHFdt(1:end - 1));
                
                % Composite
                W.primary(p).I_pPk = max(abs(W.primary(p).waveform.i_p));
                W.primary(p).I_pRMS = sqrt(W.primary(p).I_pLFRMS^2 + W.primary(p).I_pHFRMS^2);
                W.primary(p).V_pRMS = sqrt(W.primary(p).V_pLFRMS^2 + W.primary(p).V_pHFRMS^2);
                W.primary(p).VA = W.primary(p).I_pRMS*W.primary(p).V_pRMS;
                P.P_t = P.P_t + W.primary(p).VA;
                temp = [W.primary(p).waveform.i_p, W.primary(p).waveform.i_p(1)];
                W.primary(p).waveform.di_pdt = diff(temp)/T.dt;
                W.primary(p).di_pdt_RMS = sqrt(W.primary(p).di_pLFdt_RMS^2 + W.primary(p).di_pHFdt_RMS^2);
            end
        else
            % LF
            iwf = W.primary.waveform.i_pLF;
            vwf = W.primary.waveform.v_pLF;
            W.primary.I_pLFPk = max(abs(iwf));
            W.primary.I_pLFRMS = computePERMS(iwf);
            W.primary.V_pLFRMS = computePERMS(vwf);
            W.primary.VALF = W.primary.I_pLFRMS*W.primary.V_pLFRMS;
            P.P_tLF = P.P_tLF + W.primary.VALF;
            temp = [W.primary.waveform.i_pLF, W.primary.waveform.i_pLF(1)];
            W.primary.waveform.di_pLFdt = diff(temp)/T.dt_0;
            W.primary.di_pLFdt_RMS = computePERMS(W.primary.waveform.di_pLFdt(1:end - 1));

            % HF
            iwf = W.primary.waveform.i_pHF;
            vwf = W.primary.waveform.v_pHF;
            W.primary.I_pHFPk = max(abs(iwf));
            W.primary.I_pHFRMS = computePERMS(iwf);
            W.primary.V_pHFRMS = computePERMS(vwf);
            W.primary.VAHF = W.primary.I_pHFRMS*W.primary.V_pHFRMS;
            P.P_tHF = P.P_tHF + W.primary.VAHF;
            temp = [W.primary.waveform.i_pHF, W.primary.waveform.i_pHF(1)];
            W.primary.waveform.di_pHFdt = diff(temp)/T.dt_s;
            W.primary.di_pHFdt_RMS = computePERMS(W.primary.waveform.di_pHFdt(1:end - 1));

            % Composite
            W.primary.I_pPk = max(abs(W.primary.waveform.i_p));
            W.primary.I_pRMS = sqrt(W.primary.I_pLFRMS^2 + W.primary.I_pHFRMS^2);
            W.primary.V_pRMS = sqrt(W.primary.V_pLFRMS^2 + W.primary.V_pHFRMS^2);
            W.primary.VA = W.primary.I_pRMS*W.primary.V_pRMS;
            P.P_t = P.P_t + W.primary.VA;
            temp = [W.primary.waveform.i_p, W.primary.waveform.i_p(1)];
            W.primary.waveform.di_pdt = diff(temp)/T.dt;
            W.primary.di_pdt_RMS = sqrt(W.primary.di_pLFdt_RMS^2 + W.primary.di_pHFdt_RMS^2);
        end

        if nws > 1
            for s = 1:nws
                % LF
                iwf = W.secondary(s).waveform.i_sLF;
                vwf = W.secondary(s).waveform.v_sLF;
                W.secondary(s).I_sLFPk = max(abs(iwf));
                W.secondary(s).I_sLFRMS = computePERMS(iwf);
                W.secondary(s).V_sLFRMS = computePERMS(vwf);
                W.secondary(s).VALF = W.secondary(s).I_sLFRMS*W.secondary(s).V_sLFRMS;
                P.P_tLF = P.P_tLF + W.secondary(s).VALF;
                temp = [W.secondary(s).waveform.i_sLF, W.secondary(s).waveform.i_sLF(1)];
                W.secondary(s).waveform.di_sLFdt = diff(temp)/T.dt_0;
                W.secondary(s).di_sLFdt_RMS = computePERMS(W.secondary(s).waveform.di_sLFdt(1:end - 1));
                
                % HF
                iwf = W.secondary(s).waveform.i_sHF;
                vwf = W.secondary(s).waveform.v_sHF;
                W.secondary(s).I_sHFPk = max(abs(iwf));
                W.secondary(s).I_sHFRMS = computePERMS(iwf);
                W.secondary(s).V_sHFRMS = computePERMS(vwf);
                W.secondary(s).VAHF = W.secondary(s).I_sHFRMS*W.secondary(s).V_sHFRMS;
                P.P_tHF = P.P_tHF + W.secondary(s).VAHF;
                temp = [W.secondary(s).waveform.i_sHF, W.secondary(s).waveform.i_sHF(1)];
                W.secondary(s).waveform.di_sHFdt = diff(temp)/T.dt_s;
                W.secondary(s).di_sHFdt_RMS = computePERMS(W.secondary(s).waveform.di_sHFdt(1:end - 1));
                
                % Composite
                W.secondary(s).I_sPk = max(abs(W.secondary(s).waveform.i_s));
                W.secondary(s).I_sRMS = sqrt(W.secondary(s).I_sLFRMS^2 + W.secondary(s).I_sHFRMS^2);
                W.secondary(s).V_sRMS = sqrt(W.secondary(s).V_sLFRMS^2 + W.secondary(s).V_sHFRMS^2);
                W.secondary(s).VA = W.secondary(s).I_sRMS*W.secondary(s).V_sRMS;
                P.P_t = P.P_t + W.secondary(s).VA;
                temp = [W.secondary(s).waveform.i_s, W.secondary(s).waveform.i_s(1)];
                W.secondary(s).waveform.di_sdt = diff(temp)/T.dt;
                W.secondary(s).di_sdt_RMS = sqrt(W.secondary(s).di_sLFdt_RMS^2 + W.secondary(s).di_sHFdt_RMS^2);
            end
        else
            % LF
            iwf = W.secondary.waveform.i_sLF;
            vwf = W.secondary.waveform.v_sLF;
            W.secondary.I_sLFPk = max(abs(iwf));
            W.secondary.I_sLFRMS = computePERMS(iwf);
            W.secondary.V_sLFRMS = computePERMS(vwf);
            W.secondary.VALF = W.secondary.I_sLFRMS*W.secondary.V_sLFRMS;
            P.P_tLF = P.P_tLF + W.secondary.VALF;
            temp = [W.secondary.waveform.i_sLF, W.secondary.waveform.i_sLF(1)];
            W.secondary.waveform.di_sLFdt = diff(temp)/T.dt_0;
            W.secondary.di_sLFdt_RMS = computePERMS(W.secondary.waveform.di_sLFdt(1:end - 1));

            % HF
            iwf = W.secondary.waveform.i_sHF;
            vwf = W.secondary.waveform.v_sHF;
            W.secondary.I_sHFPk = max(abs(iwf));
            W.secondary.I_sHFRMS = computePERMS(iwf);
            W.secondary.V_sHFRMS = computePERMS(vwf);
            W.secondary.VAHF = W.secondary.I_sHFRMS*W.secondary.V_sHFRMS;
            P.P_tHF = P.P_tHF + W.secondary.VAHF;
            temp = [W.secondary.waveform.i_sHF, W.secondary.waveform.i_sHF(1)];
            W.secondary.waveform.di_sHFdt = diff(temp)/T.dt_s;
            W.secondary.di_sHFdt_RMS = computePERMS(W.secondary.waveform.di_sHFdt(1:end - 1));

            % Composite
            W.secondary.I_sPk = max(abs(W.secondary.waveform.i_s));
            W.secondary.I_sRMS = sqrt(W.secondary.I_sLFRMS^2 + W.secondary.I_sHFRMS^2);
            W.secondary.V_sRMS = sqrt(W.secondary.V_sLFRMS^2 + W.secondary.V_sHFRMS^2);
            W.secondary.VA = W.secondary.I_sRMS*W.secondary.V_sRMS;
            P.P_t = P.P_t + W.secondary.VA;
            temp = [W.secondary.waveform.i_s, W.secondary.waveform.i_s(1)];
            W.secondary.waveform.di_sdt = diff(temp)/T.dt;
            W.secondary.di_sdt_RMS = sqrt(W.secondary.di_sLFdt_RMS^2 + W.secondary.di_sHFdt_RMS^2);
        end
    end
end