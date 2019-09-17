%% distributeWindingWaveformsLF
% Distribute csv data amongst appropriate windings, assuming proper csv
% structure in original files:
% 
% Primary:
% Time    vp1    ip1    vp2    ip2   ...
% 
% Secondary:
% Time    vs1    is1    vs2    is2   ...
% 
% Handles LF oscillation with HF switching waveform superimposed.

function [T, Winding] = distributeWindingWaveformsLF(Winding, ptemp, stemp, tVec, T, f_s, f_0)
    thisP = Winding.primary;
    thisS = Winding.secondary;
    [~, nwp] = size(thisP);
    [~, nws] = size(thisS);
    
    time = linspace(0, 1/f_s, 1025);
    time0 = linspace(0, 1/f_0, 1025);
    timeBig = 0:time(2) - time(1):time0(end);
    
    [ptemp, t] = cleanupWaveforms(ptemp, tVec, time0(end));
    [stemp, ~] = cleanupWaveforms(stemp, tVec, time0(end));
    
    tVec = t;
    
    % separate and distribute waveforms to structures
    if nwp == 1
        vwf = ptemp(:, 2);
        iwf = ptemp(:, 3);
        [vlf, vhf] = separateWaveforms(tVec, vwf, timeBig, time, time0);
        thisP.waveform.v_p = interp1(tVec, vwf, timeBig);
        thisP.waveform.v_pLF = vlf;
        thisP.waveform.v_pHF = vhf;
        [ilf, ihf] = separateWaveforms(tVec, iwf, timeBig, time, time0);
        thisP.waveform.i_p = interp1(tVec, iwf, timeBig);
        thisP.waveform.i_pLF = ilf;
        thisP.waveform.i_pHF = ihf;
        
        % avoid numerical errors and assume periodic point for fresh interps
        thisP.waveform.v_p(end) = thisP.waveform.v_p(1);
        thisP.waveform.i_p(end) = thisP.waveform.i_p(1);
    else
        for i = 1:nwp
            offset = 2*i;
            vwf = ptemp(:, offset);
            iwf = ptemp(:, offset + 1);
            [vlf, vhf] = separateWaveforms(tVec, vwf, timeBig, time, time0);
            thisP(i).waveform.v_p = interp1(tVec, vwf, timeBig);
            thisP(i).waveform.v_pLF = vlf;
            thisP(i).waveform.v_pHF = vhf;
            [ilf, ihf] = separateWaveforms(tVec, iwf, timeBig, time, time0);
            thisP(i).waveform.i_p = interp1(tVec, iwf, timeBig);
            thisP(i).waveform.i_pLF = ilf;
            thisP(i).waveform.i_pHF = ihf;
            thisP(i).waveform.v_p(end) = thisP(i).waveform.v_p(1);
            thisP(i).waveform.i_p(end) = thisP(i).waveform.i_p(1);
        end
    end

    if nws == 1
        vwf = stemp(:, 1);
        iwf = stemp(:, 2);
        [vlf, vhf] = separateWaveforms(tVec, vwf, timeBig, time, time0);
        thisS.waveform.v_s = interp1(tVec, vwf, timeBig);
        thisS.waveform.v_sLF = vlf;
        thisS.waveform.v_sHF = vhf;
        [ilf, ihf] = separateWaveforms(tVec, iwf, timeBig, time, time0);
        thisS.waveform.i_s = interp1(tVec, iwf, timeBig);
        thisS.waveform.i_sLF = ilf;
        thisS.waveform.i_sHF = ihf;
        thisS.waveform.v_s(end) = thisS.waveform.v_s(1);
        thisS.waveform.i_s(end) = thisS.waveform.i_s(1);
    else
        for i = 1:nws
            offset = 2*i - 1;
            vwf = stemp(:, offset);
            iwf = stemp(:, offset + 1);
            [vlf, vhf] = separateWaveforms(tVec, vwf, timeBig, time, time0);
            thisS(i).waveform.v_s = interp1(tVec, vwf, timeBig);
            thisS(i).waveform.v_sLF = vlf;
            thisS(i).waveform.v_sHF = vhf;
            [ilf, ihf] = separateWaveforms(tVec, iwf, timeBig, time, time0);
            thisS(i).waveform.i_s = interp1(tVec, iwf, timeBig);
            thisS(i).waveform.i_sLF = ilf;
            thisS(i).waveform.i_sHF = ihf;
            thisS(i).waveform.v_s(end) = thisS(i).waveform.v_s(1);
            thisS(i).waveform.i_s(end) = thisS(i).waveform.i_s(1);
        end
    end
    
    T.t_s = time;
    T.dt_s = T.t_s(2) - T.t_s(1);
    T.t_0 = time0;
    T.dt_0 = T.t_0(2) - T.t_0(1);
    T.t = timeBig;
    T.dt = T.t(2) - T.t(1);

    Winding.primary = thisP;
    Winding.secondary = thisS;
    
    % separateWaveforms
    % separates low- and high-frequency portions of a waveform using a time
    % vector `tBig` with length equal to low-frequency period and sample rate
    % equal to high-frequency sample rate (after interpolation).
    function [lf, hf] = separateWaveforms(tOrig, wfOrig, tBig, t, t0)
        [~, idxTsH] = min(abs(tBig - t(end)/2));
        nHPeriods = ceil(numel(tBig)/idxTsH);
        p = interp1(tOrig, wfOrig, tBig);
        [LF, idxMax] = max(abs(p(1:idxTsH)));
        tLF = tBig(idxMax);
        
        for pt = 2:nHPeriods
            if pt < nHPeriods
                period = (pt - 1)*idxTsH + 1:pt*idxTsH;
                [maxpoint, idxMax] = max(abs(p(period)));
                LF = [LF, maxpoint];
                tLF = [tLF, tBig((pt - 1)*idxTsH + 1 + idxMax)];
            else
                period = (pt - 1)*idxTsH + 1:numel(p);
                [maxpoint, idxMax] = max(abs(p(period)));
                LF = [LF, maxpoint];
                tLF = [tLF, tBig((pt - 1)*idxTsH + 1 + idxMax)];
            end
        end
        
        % enforce periodicity and remove offset from LF waveform
        mid = (max(LF) + min(LF))/2;
        midSlope = (LF(1) - LF(end))/(tLF(1) + t0(end) - tLF(end));
        missStart = LF(1) - tLF(1)*midSlope;
        missEnd = LF(end) + (t0(end) - tLF(end))*midSlope;
        
        LF = [missStart, LF, missEnd];
        LF = LF - mid; % DC offset
        tLF = [0, tLF, t0(end)];
        
        newLF = interp1(tLF, LF, tBig);
        lf = interp1(tBig, newLF, t0);
        
        % check for complete vector (might not have gotten exact period)
        if isnan(lf(end))
            lf(end) = lf(1);
        end
        
        % assign one period of (now average peak) switching waveform
        newHF = p - sign(p).*newLF;
        hf = newHF(1:numel(t));
        
        % enforce periodicity in HF waveform by circular shift before returning
        if ~isequal(hf(1), hf(end))
            hf = circshift(hf, -1);
        end
    end
end