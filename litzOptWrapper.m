%% litzOptWrapper
% A complete wrapper for Dartmouth's LitzOpt tool.  Accepts a transformer
% structure and a time structure corresponding to it.  Returns transformer
% structure unedited for now.  May include results later.
% TODO: compute eta instead of declaring it

function T = litzOptWrapper(T, Time)
    global RHO_CU
    
    writeLitzFile(T, Time);
    
    % Run lightly edited version of litzopt
    litzoptEdit();
    
    %% Nested Functions

    % writeLitzFile
    % Create document by unpacking transformer structure.
    function writeLitzFile(T, Time)
        [C, P, W] = unpackTransformer(T);
        
        % numbers of windings
        nw = P.N_w;
        nwp = P.N_wp;
        nws = P.N_ws;
        
        % per-winding parameters with assumptions:
        % * Full bobbin
        % * Even distribution of window area for all windings
        % Note: centers have x and y reversed from directions! (e.g x: h, y: b)
        Wgap = W.t_ins;
        heightVal = (C.bobbin.height - Wgap*(nw - 1))/nw;
        breadthVal = C.bobbin.breadth;
        bht = C.window.height - C.bobbin.height; % offset, h direction
        bbt = (C.window.breadth - C.bobbin.breadth)/2; % offset, b direction
        hBase = bht;
        rLegC = (C.d_center + C.d_center2)/4; % circular center leg radius
        Istr = {};
        
        % packing factor considering bobbin and litz wire porosity factor
        eta_litz = 0.75;
        eta_bobbin = (C.bobbin.breadth*C.bobbin.height)/(C.window.breadth*C.window.height);
        eta = eta_litz*eta_bobbin;
        
        Nstr = 'N = [';
        hstr = 'a = [';
        wstr = 'b = [';
        cxstr = 'xzero = [';
        cystr = 'yzero = [';
        lenstr = 'len = [';
        endstr = '];\n';
        w = 0;
        
        if nwp > 1
            for p = 1:nwp                
                Nstr = [Nstr, num2str(W.primary(p).N, '%i'), ', '];
                hstr = [hstr, num2str(heightVal, '%g'), ', '];
                wstr = [wstr, num2str(breadthVal, '%g'), ', '];
                cxstr = [cxstr, num2str(hBase + heightVal/2, '%g'), ', '];
                cystr = [cystr, num2str(bbt + breadthVal/2, '%g'), ', '];
                lenstr = [lenstr, ...
                          num2str(2*pi*(rLegC + hBase + heightVal/2), '%g'), ...
                          ', '];
                Istr{w + 1, 1} = [sprintf('I(%d, :) = [', w + 1), stringifyVector(W.primary(p).waveform.i_p9(1:end - 1)), endstr];
                
                hBase = hBase + heightVal + Wgap;
                w = w + 1;
            end
        else
            Nstr = [Nstr, num2str(W.primary.N, '%i'), ', '];
            hstr = [hstr, num2str(heightVal, '%g'), ', '];
            wstr = [wstr, num2str(breadthVal, '%g'), ', '];
            cxstr = [cxstr, num2str(hBase + heightVal/2, '%g'), ', '];
            cystr = [cystr, num2str(bbt + breadthVal/2, '%g'), ', '];
            lenstr = [lenstr, ...
                      num2str(2*pi*(rLegC + hBase + heightVal/2), '%g'), ...
                      ', '];
            Istr{w + 1, 1} = [sprintf('I(%d, :) = [', w + 1), stringifyVector(W.primary.waveform.i_p9(1:end - 1)), endstr];
            
            hBase = hBase + heightVal + Wgap;
            w = w + 1;
        end

        if nws > 1
            for s = 1:nws
                Nstr = [Nstr, num2str(W.secondary(s).N, '%i')];
                hstr = [hstr, num2str(heightVal, '%g')];
                wstr = [wstr, num2str(breadthVal, '%g')];
                cxstr = [cxstr, num2str(hBase + heightVal/2, '%g')];
                cystr = [cystr, num2str(bbt + breadthVal/2, '%g')];
                lenstr = [lenstr, ...
                          num2str(2*pi*(rLegC + hBase + heightVal/2), '%g')];
                Istr{w + 1, 1} = [sprintf('I(%d, :) = [', w + 1), stringifyVector(W.secondary(s).waveform.i_s9(1:end - 1)), endstr];
                
                if w < nw
                    Nstr = [Nstr, ', '];
                    hstr = [hstr, ', '];
                    wstr = [wstr, ', '];
                    cxstr = [cxstr, ', '];
                    cystr = [cystr, ', '];
                    lenstr = [lenstr, ', '];
                    
                    hBase = hBase + heightVal + Wgap;
                    w = w + 1;
                else
                    Nstr = [Nstr, endstr];
                    hstr = [hstr, endstr];
                    wstr = [wstr, endstr];
                    cxstr = [cxstr, endstr];
                    cystr = [cystr, endstr];
                    lenstr = [lenstr, endstr];
                end
            end
        else
            Nstr = [Nstr, num2str(W.secondary.N, '%i'), endstr];
            hstr = [hstr, num2str(heightVal, '%g'), endstr];
            wstr = [wstr, num2str(breadthVal, '%g'), endstr];
            cxstr = [cxstr, num2str(hBase + heightVal/2, '%g'), endstr];
            cystr = [cystr, num2str(bbt + breadthVal/2, '%g'), endstr];
            lenstr = [lenstr, ...
                      num2str(2*pi*(rLegC + hBase + heightVal/2), '%g'), ...
                      endstr];
            Istr{w + 1, 1} = [sprintf('I(%d, :) = [', w + 1), stringifyVector(W.secondary.waveform.i_s9(1:end - 1)), endstr];
        end
        
        [fileName, pathName] = uiputfile([T.core.name, '.m']);
        fid = fopen(fullfile(pathName, fileName), 'w');
        fprintf(fid, 'rhoc = %g;\n', RHO_CU);
        fprintf(fid, 'awg = 32:2:48;\n');
        fprintf(fid, 'Ximage = 5;\n');
        fprintf(fid, 'Yimage = 5;\n');
        fprintf(fid, 'xdiv = 20;\n');
        fprintf(fid, 'ydiv = 20;\n');
        fprintf(fid, 'h = %g;\n', C.window.height);
        fprintf(fid, 'bw = %g;\n', C.window.breadth);
        fprintf(fid, 'hb = %g;\n', C.bobbin.height);
        fprintf(fid, 'bb = %g;\n', C.bobbin.breadth);
        fprintf(fid, 'gaplen = %g;\n', 0); % maybe change this
        fprintf(fid, 'insulbuild = ''s'';\n'); % ... and this
        fprintf(fid, 'gaploc = ''No Gap'';\n'); % ... and this
        fprintf(fid, 'userfp = %g;\n', eta);
        fprintf(fid, Nstr);
        fprintf(fid, hstr);
        fprintf(fid, wstr);
        fprintf(fid, cxstr);
        fprintf(fid, cystr);
        fprintf(fid, lenstr);
        fprintf(fid, ['dt = [', stringifyVector(diff(Time.t9)), endstr]);
        
        % original design parameters (not evaluated)
        fprintf(fid, 'orig_wiresize = [0, 0];\n');
        fprintf(fid, 'orig_numstrands = [0, 0];\n');
        
        % current vectors
        for idx = 1:size(Istr, 1)
            fprintf(fid, Istr{idx});
        end
    end

    % unpackTransformer
    % Unpacks core, properties, and winding for ease of access (and brevity).
    function [C, P, W] = unpackTransformer(T)
        W = T.winding;
        C = T.core;
        P = T.properties;
    end

    % stringifyVector
    % Adds elements of vector to string with comma separation
    function s = stringifyVector(v)
        s = '';
        
        for idx = 1:length(v) - 1
            s = [s, num2str(v(idx), '%g'), ', '];
        end
        
        s = [s, num2str(v(end), '%g')];
    end
end