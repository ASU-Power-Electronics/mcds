%% distributeWindingWaveforms
% Distribute csv data amongst appropriate windings, assuming proper csv
% structure in original files:
% 
% Primary:
% |Time    vp1    ip1    vp2    ip2   ...|
% 
% Secondary:
% |Time    vs1    is1    vs2    is2   ...|

function Winding = distributeWindingWaveforms(Winding, ptemp, stemp, tVec, time)
    thisP = Winding.primary;
    thisS = Winding.secondary;
    [~, nwp] = size(thisP);
    [~, nws] = size(thisS);
    Ts = time(end);
    
    [ptemp, t] = cleanupWaveforms(ptemp, tVec, Ts);
    [stemp, ~] = cleanupWaveforms(stemp, tVec, Ts);
    
    tVec = t;

    if nwp == 1
        thisP.waveform.v_p = interp1(tVec, ptemp(:, 2), time);
        thisP.waveform.i_p = interp1(tVec, ptemp(:, 3), time);
        % avoid numerical errors and assume periodic point
        thisP.waveform.v_p(end) = thisP.waveform.v_p(1);
        thisP.waveform.i_p(end) = thisP.waveform.i_p(1);
    else
        for i = 1:nwp
            offset = 2*i;
            thisP(i).waveform.v_p = interp1(tVec, ptemp(:, offset), time);
            thisP(i).waveform.i_p = interp1(tVec, ptemp(:, offset + 1), time);
            thisP(i).waveform.v_p(end) = thisP(i).waveform.v_p(1);
            thisP(i).waveform.i_p(end) = thisP(i).waveform.i_p(1);
        end
    end

    if nws == 1
        thisS.waveform.v_s = interp1(tVec, stemp(:, 1), time);
        thisS.waveform.i_s = interp1(tVec, stemp(:, 2), time);
        thisS.waveform.v_s(end) = thisS.waveform.v_s(1);
        thisS.waveform.i_s(end) = thisS.waveform.i_s(1);
    else
        for i = 1:nws
            offset = 2*i - 1;
            thisS(i).waveform.v_s = interp1(tVec, stemp(:, offset), time);
            thisS(i).waveform.i_s = interp1(tVec, stemp(:, offset + 1), time);
            thisS(i).waveform.v_s(end) = thisS(i).waveform.v_s(1);
            thisS(i).waveform.i_s(end) = thisS(i).waveform.i_s(1);
        end
    end

Winding.primary = thisP;
Winding.secondary = thisS;
end