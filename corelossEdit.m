%% CORELOSS Caclulate power loss in a ferrite core.
%   
%   CORELOSS(time, B, alpha, beta, k, suppressFlag)
%
%   time and B are a piecewise linear series of data points
%   time = vector of successive time values corresponding to flux values, in seconds.
%   B = vector of flux densities corresponding to time vector values, in tesla.
%   alpha, beta, k = Steinmetz parameters.  Must be for use with MKS units
%   suppressFlag: 0 to display input plot, loss  and error information in the
%                   command window.
%                 1 to suppress all output and only return the loss value.
%                   If an error occurs, CORELOSS will return -1.
%             
%   CORELOSS will calculate the power loss in your specified ferrite core.
%   The loss value is in W/m^3.
%
% Created by Charles Sullivan, Kapil Venkatachalam, Jens Czogalla
% at Thayer School of Engineering, Dartmouth College.
% Modified by Kip Benson 3/03 (error checking on input data and improved
% command line interface).
% Corrected 11/05 C. Sullivan:  fixed minor loop handling bug ("minortime" was mistyped as minorime")
%
% Updated 9/19 R. S. Mongrain:  Modernized for speed and memory usage.  Comments
% cleaned up, variables given self-documenting names.  Full implementation of
% iGSE used instead of PWL approximation.

function y = corelossEdit(tvs, B, alpha, beta, k, suppressFlag)
    % tvs is what was called "time" in help info above.
    % check data for errors
	error = 'None';
	
	% makes sure times are successive
	if any(diff(tvs) < 0)
    	error = 'Time data points must be successive.';
	end
    
    % make sure vectors are same length
    if ~isequal(length(tvs), length(B))
        error = 'Time and flux vectors must have same length';
    end
    
    if ~isequal(B(1), B(end))
        error = 'Since the PWL input is periodic, the first flux value must equal the last';
    end
    
    % R. Scott Mongrain - force row vectors to avoid indexing errors
    if size(tvs, 1) > size(tvs, 2)
        tvs = tvs';
    end
    
    if size(B, 1) > size(B, 2)
        B = B';
    end

    % calculate core loss if there is no error with the input data
    if suppressFlag == 0

        if strcmp(error, 'None')
            Pcore = 1000*gsepwl(tvs', B, alpha, beta, k);
            sprintf('Core Loss:  %g [W/m^3]\n', Pcore)
            y = Pcore;

            % plot the user's input
            pwlPlot = figure;
            pwlPlot.Visible = 'on';
            plot(tvs, B, 'b.-');
            xlabel('Time (s)');
            ylabel('Flux Density (T)');
            title('Plot of Input Data');
        else
            sprintf('Error:  %s\n', error)
            y = -1;
        end

    elseif suppressFlag == 1	% suppress output
        if strcmp(error, 'None')
            y = 1000*gsepwl(tvs', B, alpha, beta, k);
        else
            y = -1;	% error occured
        end

    end


    % gsepwl
    % to calculate core loss per unit volume using iGSE
    function p = gsepwl(t, B, alpha, beta, k)
                    % a is fraction of Bpp used, 1 - a is fraction of original
        a = 1.0;    % a=1 for iGSE
                    % a=0 for GSE
        T = t(end) - t(1);	% total time of PWL period
%         ki = k/((2^(beta + 1))*pi^(alpha - 1)*(0.2761 + 1.7061/(alpha + 1.354)));
        ki = k/(2^(beta - 1)*pi^(alpha - 1)*integral(@(theta) abs(cos(theta)).^alpha, 0, 2*pi));
        [B, t, Bs, ts] = splitLoop(B, t);	% split waveform into major and minor loops
        pseg(1) = calcSeg(t, B, alpha, beta, ki, a);
        dt(1) = t(end) - t(1);

        for seg = 1:numel(ts)
            pseg(seg + 1) = calcSeg(ts{seg}, Bs{seg}, alpha, beta, ki, a);
            tseg = ts{seg};
            dt(seg + 1) = tseg(end) - tseg(1);
        end

        p = sum(pseg.*dt)/T;
    end

    % splitLoop
    % Loop Splitter by Kapil Venkatachalam
    % Modified on June 2nd, 2002
    %
    % This function takes a piecewise linear current waveform corresponding to a
    % B-H loop and splits it into smaller waveforms which have the same starting
    % and ending value.
    %
    % Inputs:
    %
    % `B`:  B vector corresponding to the B-H loop.
    % `t`:  Time vector corresponding to the B-H loop. 
    %
    % Outputs:
    %
    % `majorloop`:  B values corresponding to the major B-H loop.
    % `majortime`:  time vector corresponding to the major B-H loop.
    % `minorloop`:  B values corresponding to the different minor B-H loops.
    % `minortime`:  time vector corresponding to the different minor B-H loops.

    function [majorloop, majortime, minorloop, minortime] = splitLoop(B, t)
        % Reshaping the waveforms and identifying the peak point
        % Identifies the lowest point in `B`.       
        [Bmod(1), idxMin] = min(B); % adds the lowest point as the first value.

        % `Bmod` is a vector storing the shifted values of `B`.
        % `tShift` is a vector storing the cumulative shifted values of `t`.        
        tShift = 0;         % adds the coresponding time value the first value.
        tDiff = diff(t);	% adjusts for the times corresponding to the values.

        % Stores all value of `B` (from the min to the end) in `Bmod`.
        Bmod = [Bmod, B(idxMin + 1:end)];
        tShift = [tShift; cumsum(tDiff(idxMin:end))];

        % Stores all value of `B` (from the beginning to the min) in `Bmod`.
        Bmod = [Bmod, B(2:idxMin)];
        tShift = [tShift; cumsum(tDiff(1:idxMin - 1))];

        % Finds the position of the peak value
        % Constant defined to identify the position of the peak value.
        [~, idxMax] = max(Bmod);

        % Defining variables to keep track of values in the vectors
        idxBmod = 2;	% current `Bmod` index
        idxMinorLoop = 1;  % current `minorLoopVals` index

        minorLoopVals = {};     % Defining cells for minorloop.
        minorLoopTimes = {};	% Defining cells for minortime.

        % `majorVals` stores the majorloop values from `Bmod`.
        % `majorTimes` stores the corresponding time values from `tShift`.
        majorVals(1) = Bmod(1);
        majorTimes(1) = tShift(1);

        % Splits the rising portion of the waveform
        while idxBmod <= idxMax
            % Check to see if the waveform is rising
            isRising = Bmod(idxBmod) >= Bmod(idxBmod - 1);
            
            if isRising
                majorVals = [majorVals, Bmod(idxBmod)];
                majorTimes = [majorTimes, tShift(idxBmod)];
                idxBmod = idxBmod + 1;
            else
                % Check for minor loop in the rising part
                if numel(minorLoopVals) < idxMinorLoop
                    minorLoopVals{idxMinorLoop, 1} = Bmod(idxBmod - 1);
                    minorLoopTimes{idxMinorLoop, 1} = tShift(idxBmod - 1);
                else
                    minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod - 1)];
                    minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod - 1)];
                end

                % R. Scott Mongrain - vectorization
                idxMajorMax = find(Bmod(idxBmod:end) >= max(majorVals), 1) + idxBmod - 1;
                minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod:idxMajorMax - 1)];
                minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod:idxMajorMax - 1)'];
                idxBmod = idxMajorMax;
                
%                 % Repeating the process till the minor loop ends
%                 while Bmod(idxBmod) < max(majorVals)
%                     minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod)];
%                     minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod)];
%                     idxBmod = idxBmod + 1;
%                 end

                % Calculating the slope of the rising edge of the minor loop.
                slope = (Bmod(idxBmod - 1) - Bmod(idxBmod))/(tShift(idxBmod - 1) - tShift(idxBmod)); 

                % Makes the last element of the minor loop same as the maximum value of `majorVals`.
                minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, max(majorVals)];

                % Computes the value of time of the point which was stored last in `minorLoopVals`.
                stemp = (max(majorVals) - Bmod(idxBmod - 1))/slope + tShift(idxBmod - 1);	% Calculating the time value.
                minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, stemp];
                majorVals = [majorVals, max(majorVals)];	% The last point in `majorVals` is repeated for continuity.
                majorTimes = [majorTimes, stemp];           % The time value is also repeated.
                idxMinorLoop = idxMinorLoop + 1;
            end  
        end

        % `majorVals` now stores the rising part of the major loop.
        % `majorTimes` stores the corresponding time values of the major loop.
        % `minorLoopVals` stores the minor loops in the rising part of the waveform.
        % `minorLoopTimes` stores the corresponding time values of the minor loops.

        % Splits the falling portion of the waveform        
        while idxBmod <= length(Bmod)
            isFalling = Bmod(idxBmod) <= Bmod(idxBmod - 1);
            
            if isFalling	% Compares adjacent values of `Bmod`.
                majorVals = [majorVals, Bmod(idxBmod)];
                majorTimes = [majorTimes, tShift(idxBmod)];
                idxBmod = idxBmod + 1;
            else
                %Check for minor loop in the falling part.
                temp = Bmod(idxBmod - 1);	% Temporary variable to store last value in `Bmod`.
                
                if numel(minorLoopVals) < idxMinorLoop
                    minorLoopVals{idxMinorLoop, 1} = Bmod(idxBmod-1);
                    minorLoopTimes{idxMinorLoop, 1} = tShift(idxBmod-1);
                else
                    minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod-1)]; 
                    minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod-1)];
                end
                
                isMinLoop = 1;  % Binary flag indicating minor loop, added by Jens Czogalla

                idxLastGT = find(Bmod(idxBmod:end) <= temp, 1) + idxBmod - 1;
                idxEndPoint = min(length(Bmod), idxLastGT);
                minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod:idxEndPoint)];
                minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod:idxEndPoint)'];
                idxBmod = idxEndPoint;
                
%                 while idxBmod <= length(Bmod) && Bmod(idxBmod) > temp
%                     minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, Bmod(idxBmod)];
%                     minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, tShift(idxBmod)];
%                     isMinLoop = 1;
%                     idxBmod = idxBmod + 1;
%                 end

                % Repeating the process till the minor loop ends
                % by Jens Czogalla
                while idxBmod <= length(Bmod) && isMinLoop
                    % Calculating the slope of the rising edge of the minor loop.
                    slope = (Bmod(idxBmod - 1) - Bmod(idxBmod))/(tShift(idxBmod - 1) - tShift(idxBmod));

                    % Makes the first element of the minor loop same as the last element in `majorVals`.
                    minorLoopVals{idxMinorLoop, 1} = [minorLoopVals{idxMinorLoop, 1}, temp];

                    % Computes the value of time of the point which was stored last in `minorLoopVals`.
                    qtemp = ((temp - Bmod(idxBmod - 1)) / slope) + tShift(idxBmod - 1);	% Calculating the time value.
                    minorLoopTimes{idxMinorLoop, 1} = [minorLoopTimes{idxMinorLoop, 1}, qtemp];

                    if ~isequal(Bmod(idxBmod), temp)
                        majorVals = [majorVals, temp];      % The last point in `majorVals` is repeated for continuity.
                        majorTimes = [majorTimes, qtemp];	% The time value is also repeated.
                    end

                    isMinLoop = 0;
                    idxMinorLoop = idxMinorLoop + 1;
                end
            end   
        end
        
        % Removal of repetition of points in `majorVals` and adjusting the values
        majorloop = majorVals(1);                 
        majortime = majorTimes(1);    
        tAdj = 0;	% Variable that adjusts time after the point of repetition.

        % by Jens Czogalla
        for idx = 2:length(majorVals)
            if ~isequal(majorVals(idx - 1), majorVals(idx))	% Checking for repeated points.
            	majortime = [majortime, majorTimes(idx) - tAdj];
            	majorloop = [majorloop, majorVals(idx)];    
            else
                tAdj = tAdj + majorTimes(idx) - majorTimes(idx - 1);	% Adjusts the time to compensate for repetitions.
            end   
        end

        % Finding the number of minor loops to be used for checking sub loops later
        nMinLoops = numel(minorLoopVals);   % Variable to keep track of the number of minor loops.

        % Recursion
        % Checks if any portion of the  split waveform has subloops.
        % If so repeats the above process to make it a single loop.
        minorloop = {};
        minortime = {};

        for curLoop = 1:nMinLoops
            curMinLoop = minorLoopVals{curLoop};	% flux
            curMinTimes = minorLoopTimes{curLoop};	% time
            
            % R. Scott Mongrain - force column vector for consistency
            if size(curMinTimes, 2) > size(curMinTimes, 1)
                curMinTimes = curMinTimes';
            end

            if isMinorLoop(curMinLoop, curMinTimes)
                [majLoop, majTime, minLoop, minTime] = splitLoop(curMinLoop, curMinTimes);    
                if ~isempty(majLoop)
                    minorloop{length(minorloop) + 1} = majLoop;
                    minortime{length(minortime) + 1} = majTime;
                end
                
                if ~isempty(minLoop)
                    minorloop{length(minorloop) + 1} = minLoop{1};
                    minortime{length(minortime) + 1} = minTime{1};                    
                end
            else
                if ~isempty(curMinLoop)
                    minorloop{length(minorloop) + 1} = curMinLoop;
                    minortime{length(minortime) + 1} = curMinTimes;
                end
            end

        end
    end


    % calcSeg
    % calculate loss for loop segment using improved equation
    % a is fraction of Bpp used, 1 - a is fraction of original
    function pseg = calcSeg(t, B, alpha, beta, k1, a)                                    
        bma1a = (beta - alpha)*(1 - a);	% exponent of B(t)
        bmaa = (beta - alpha)*a;        % exponent of Bpp
        Bpp = max(B) - min(B);
        
        % force row vectors
        if size(B, 1) > size(B, 2)
            B = B';
        end
        
        if size(t, 1) > size(t, 2)
            t = t';
        end

        if sum(B < 0) > 0
            [t, B] = makePositive(t, B);
        end

        T = t(end) - t(1);
        deltaB = abs(diff(B));
        deltat = abs(diff(t));
        dBdt = deltaB./deltat;

        % Correct infinite slope errors - R. Scott Mongrain 08SEP19
        % Removes point singularities by taking the mean of surrounding points,
        % or by enforcing continuity of periodic signal at endpoints.
        if sum(isinf(dBdt)) > 0
            for idx = 1:length(dBdt)
                if isinf(dBdt(idx))
                    if idx > 1 && idx < length(dBdt)
                        dBdt(idx) = abs(mean([dBdt(idx - 1), dBdt(idx + 1)]));
                    else
                        if isequal(idx, 1)
                            dBdt(idx) = abs(dBdt(end));
                        else
                            dBdt(idx) = abs(dBdt(1));
                        end
                    end
                end
            end
        end

        m1 = dBdt.^(alpha - 1);
        m2 = abs((B(2:end)).^(bma1a + 1) - (B(1:end - 1)).^(bma1a + 1));
        pseg = k1/T/(bma1a + 1)*sum(m1.*m2)*Bpp^bmaa;
    end


    % makePositive
    % Convert a piecewise-linear waveform w (which must be a vector)
    % represented by the points w at times t, to a piecewise-linear waveform with 
    % points at any zero crossings.  Thus, any segment of the returned waveform is all 
    % positive or all negative.  If w is a matrix, use makepositiveM.
    function [t, w] = makePositive(t, w)
        len = length(t);
        crossing = w(1:end - 1).*w(2:end) < 0;
        loopnumber = len - 1 + sum(crossing);
        
        for idx = 1:loopnumber
           if crossing(idx)
              len = length(w);
              tcross = w(idx)*(t(idx + 1) - t(idx))/(w(idx) - w(idx + 1)) + t(idx);
              t = [t(1:idx), tcross, t(idx + 1:end)];
              w = [w(1:idx), 0, w(idx + 1:end)];
              crossing = [crossing(1:idx), 0, crossing(idx + 1:len - 1)];
           end
        end
        
        w = abs(w);
    end



    % isMinorLoop
    % Boolean-valued function to determine if waveform s over time p contains a
    % minor loop.
    function ml = isMinorLoop(wf, t)
    peak = -1;
    prevslope = 0;

        for idx = 2:length(wf)
            slope = (wf(idx - 1) - wf(idx))/(t(idx - 1) - t(idx));

            if (slope <= 0 && prevslope >= 0)
                peak = peak + 1;
            end

            if (prevslope <= 0 && slope >= 0)
                peak = peak + 1;
            end

            prevslope = slope;
        end

        if (peak > 2)
            ml = 1;
        else
            ml = 0;
        end
    end
end