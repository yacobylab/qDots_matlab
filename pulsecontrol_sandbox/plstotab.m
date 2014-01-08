function pulse = plstotab(pulse)
% pulse = plstotab(pulse)
%
% Convert 'elem' to 'tab' pulse format.
%
% Integer pulse is taken from the database.
% Fields of pulse:
% name, xval, and taurc remain unchanged.
%
% format: must be 'elem' or 'tab'. Nothing is done for 'tab';
%
% data: struct array with fields type, time, val corresponding to pulse elements
%           to be concatenated.  type is a string specifying the type of element, 
%           val and time specicy pulse voltages and times, respectively.
%           Their meaning and the format of val depend on type.
% 
% Possible type strings and corresponding interpretation of val:
% raw: insert [time; val] into pulse table.
% mark: add time' to marktab
% fill: stretch this element to make the total pulse duration equal to time.
%       Idea for future development: allow several fills, each spreading a subset.
%       Would need a second element to flush previous fill, could be fill without time.
% wait: stay at val (row vector, one entry for each channel) for duration time.
%   If val has 3 entries, third is a scaling factor for the first two.
% reload: relaod pulse at val (row vector, one entry for each channel).
%         time: [ramp time, wait time at load point, wait time at (0, 0) after load] 
% meas: measurement stage at [0, 0] for time(1), RF marker delayed by time(2) and
%       off time(3) before end of the stage.  [time(2) is lead delay,
%       time(3) is negative tail delay. 
%       Optional val(1) is the readout tag. If it is given and not nan, time 4 and 5 set its delays
%       with the sae convention as for the marker.
%       Optional val(2,3) moves the measurement point away from 0,0.  Makes
%       meas_o obsolete.
% meas_o: as meas, but measure at current voltages, not 0,0
% ramp: ramp to val (row vector, one entry for each channel) in time.  opt val(3) is multiplier
% comp: measurement compensation at val(1:2) (one for each channel) for duration time(1). 
%       Ramps voltage to target and back over time(2) and time(3) at the beginning and 
%       end of the stage, respectively. If length(val)>=4, val(3:4) are used as final value.
%       The compensation value could be determined automatically, but this feature is not 
%       implemented yet.
% adprep: adiabatic ramp along second diagonal (epsilon) from val(1) to val(2), ramp duration time.
% adread: same, going the other way.
% RFburst: make awg output sin wave: time is duration (in us).
%       val(1) = freq in MHz, val(2:3) = amplitude of x and y, 
%       optional val(4:5) = offset to 2 and 3
%       optional val(6:7) = phase for each channel
% markerburst: is a wait, with markers also firing
%       time(1)=wait/marker time, val(1:3) where to wait, (just like exch)
%       time(2:3) = fire marker this time early, fire marker this long after exch
%       val(4:7)= which marker to fire, put in 0 or 1 for C1M1,C1M2,C2M1,C2M2

% MARKTAB:
%  each column gives a marker pulse
%  [ start_time; ch1_mk1_wid ; ch1_mk2_wid ; ch_2_mk1_wid ...]
% ie,
%  [ 2 ; 0 ; 1 ; 0 ; 1 ] fires markers ch1mk1 and ch2mk2 for 1 us starting at 2us.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.



pulse = plsdefault(pulse);

% read from database, assumed to have valid format
% if ~isstruct(pulse)
%     pulse = plsdata(pulse);
% end
% while strcmp(pulse.format, 'ind')
%     pulse = plsdata.pulses(pulse.data);
% end
        
dt=-1e-9;  % Shortest meaninful length

switch pulse.format
    case 'tab'
        return;
       
    case 'elem'
        
        pulsetab = zeros(3, 0);
        marktab =  zeros(5, 0);
        comppos = [];
        fillpos = [];
        readout = [];
        readpos = [];
        
        pulsedef = pulse.data;
        
        for i = 1:length(pulsedef)

            switch pulsedef(i).type

                case 'raw'
                    pulsetab = [pulsetab, [pulsedef(i).time; pulsedef(i).val]];

                case 'mark'
                    marktab = [marktab, pulsedef(i).time'];

                case 'fill'
                    fillpos = size(pulsetab, 2);
                    filltime = pulsedef(i).time(1);
                    fillmarkpos = size(marktab,2);                    
                case 'wait'
                    if pulsedef(i).time(1) > 1e-11
                        fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        if length(pulsedef(i).val) > 2
                          pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(3)*pulsedef(i).val(1:2)', 1, 2);
                        else
                          pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(1:2)', 1, 2);
                        end
                    end

                case 'reload'
                    % If we're filling the load, push the fillpos 1 forward
                    % so we stretch the wait at the loadpos, not the ramp
                    % to the loadpos                    
                    % Ignore zero length loads
                    if pulsedef(i).time(2) > 1e-11
                      fillload = (fillpos == size(pulsetab,2));                        
                      pulsetab(1, end+(1:4)) = pulsetab(1, end) + cumsum(pulsedef(i).time([1 2 1 3]));
                      pulsetab(2:3, end+(-3:0)) = [repmat(pulsedef(i).val(1:2)', 1, 2), zeros(2)];
                      fillpos = fillpos + fillload;                    
                    end
                case 'meas_o' % offset measurement
                    pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                    pulsetab(2:3, end-1) = pulsetab(2:3,end-2);
                    pulsetab(2:3, end) = pulsetab(2:3,end-2);
                    marktab(:, end+1) = [pulsetab(1, end-2)+pulsedef(i).time(2); 0; 0; 0; pulsedef(i).time(1:3)*[1; -1; -1]];
                    if ~isempty(pulsedef(i).val)
                        readout(end+1, :) = [pulsedef(i).val, pulsetab(1, end-2) + pulsedef(i).time(4), pulsedef(i).time([1 4 5])*[1; -1; -1]];
                        readpos(end+1) = size(pulsetab, 2)-2;
                    end
                    
                case 'meas'
                    if length(pulsedef(i).val) == 3
                        mpnt = pulsedef(i).val(2:3);
                    else
                        mpnt = [0,0];
                    end
                    pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                    pulsetab(2:3, end-1) = mpnt;
                    pulsetab(2:3, end) = mpnt;
                    marktab(:, end+1) = [pulsetab(1, end-2)+pulsedef(i).time(2); 0; 0; 0; pulsedef(i).time(1:3)*[1; -1; -1]];
                    if length(pulsedef(i).val) > 0 && ~isnan(pulsedef(i).val(1))
                        readout(end+1, :) = [pulsedef(i).val(1), pulsetab(1, end-2) + pulsedef(i).time(4), pulsedef(i).time([1 4 5])*[1; -1; -1]];
                        readpos(end+1) = size(pulsetab, 2)-2;
                    end

                case 'ramp'
                    %allow for multiplies in ramps- helps get direction
                    %right
                    if length(pulsedef(i).val) ==3
                        mult = pulsedef(i).val(3);
                    else
                        mult = 1;
                    end
                    pulsetab(1, end+1) = pulsetab(1, end) + pulsedef(i).time(1);
                    pulsetab(2:3, end) = mult*pulsedef(i).val(1:2);

                case 'comp'
                    comppos = size(pulsetab, 2)+1;
                    compval  = pulsedef(i).val(1:2);

                    pulsetab(1, end+(1:4)) = pulsetab(1, end) + [0 pulsedef(i).time(2), pulsedef(i).time(1)-sum(pulsedef(i).time(2:3)), ...
                        pulsedef(i).time(1)];
                    pulsetab(2:3, end+(-3:0)) = 0;
                    if length(pulsedef(i).val) >= 4
                        pulsetab(2:3, end) = pulsedef(i).val(3:4);
                    end

                case 'adprep'
                    if pulsedef(i).time(1) > 1e-11
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)];
                        if(length(pulsedef(i).val) <= 2)
                            dir=[-1 1];
                        else
                            dir = pulsedef(i).val(3:4);
                        end
                        pulsetab(2:3, end-1) = pulsedef(i).val(1)  * dir;
                        pulsetab(2:3, end) = pulsedef(i).val(2) * dir;
                    end
                case 'adread'
                    if pulsedef(i).time(1) > 1e-11
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)];
                        if(length(pulsedef(i).val) <= 2)
                            dir=[-1 1];
                        else
                            dir = pulsedef(i).val(3:4);
                        end
                        pulsetab(2:3, end-1) = pulsedef(i).val(2)  * dir;
                        pulsetab(2:3, end) = pulsedef(i).val(1)  * dir;
                    end
                case 'RFburst' %added 2013/01/22
                    a = pulsedef(i).val;
                    a(isnan(a))=0;
                    a(1+length(a):7)=0;
                    offst = a(4:5);
                    ph = a(6:7);
%                     if length(pulsedef(i).val)>5 && all(~isnan(pulsedef(i).val(4:6)))
%                         offst = pulsedef(i).val(4:5);
%                         ph = pulsedef(i).val(6);
%                     else
%                        offst = [0 0];
%                        ph = 0;
%                     end
                    ts = 0:abs(1e6*dt):pulsedef(i).time(1);
                    if ~exist('firstburst','var') % the first burst will determine the phase for the rest of the pulse
                       firstburst = pulsetab(1,end); 
                    end
                    if length(ts) > 1 % number of time points
                        if pulsedef(i).val(1) > 1/abs(dt) % more than 1GHz
                            warning('asking for frequency above 1GHZ');
                        end
                        t_off = pulsetab(1,end)-firstburst;
                        pulsetab(1,end+(1:length(ts)))=pulsetab(1,end)+ts;
                        pulsetab(2,end+(-length(ts)+1:0)) = offst(1)+pulsedef(i).val(2)*sin(ph(1)+2*pi*pulsedef(i).val(1)*(ts+t_off));
                        pulsetab(3,end+(-length(ts)+1:0)) = offst(2)+pulsedef(i).val(3)*sin(ph(2)+2*pi*pulsedef(i).val(1)*(ts+t_off));
                    else
                        warning('asking for only one point of oscillations');
                    end
                case 'markerburst'
                    ctime = pulsetab(1,end); %current time
                    %first treat it as a wait
                    if pulsedef(i).time(1) > 1e-11
                        fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        if length(pulsedef(i).val) > 2
                            pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(3)*pulsedef(i).val(1:2)', 1, 2);
                        else
                            pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(1:2)', 1, 2);
                        end
                        % now deal with the markers
                        %marktab = [marktab, pulsedef(i).time'];
                        if length(pulsedef(i).time)>1 && all(~isnan(pulsedef(i).time(2:3)))
                           stdly = pulsedef(i).time(2); %start delay
                           edly = pulsedef(i).time(3); %end delay
                        else
                            stdly = 0; edly = 0;
                        end
                        mk = pulsedef(i).val(4:7); mk(isnan(mk))=0;
                        marktab = [marktab, [ctime-stdly, (pulsedef(i).time(1)+edly+stdly)*mk]'];
                    end
                otherwise
                    error('Invalid pulse element %i: %s.\n', i, pulsedef(i).type)
            end
        end

        if ~isempty(comppos)
            pulsetab(2:3, comppos+(1:2)) = 2*repmat(compval(1:2)', 1, 2);
        end

        %pulsetab(2:3, :) = pulsetab(2:3, :)./pinf.scale;
        %pinf = rmfield(pinf, 'scale');

        if ~isempty(fillpos)
            filltime = filltime - pulsetab(1, end);
            if filltime < 0
                pulsetab
                error('Pulse too long by %g (target %g).',-filltime,filltime+pulsetab(1,end));
            end
            pulsetab(1, fillpos+1:end) = pulsetab(1, fillpos+1:end) + filltime;
            if ~isempty(readpos)
                readout(readpos > fillpos, 2) = readout(readpos > fillpos, 2) + filltime;
            end
            marktab(1, fillmarkpos+1:end) = marktab(1, fillmarkpos+1:end) + filltime;
        end

        mask = all(abs(diff(pulsetab(2:3, :), [], 2)) < 1e-14);
        pulsetab(:, [false, mask(2:end)&mask(1:end-1)]) = [];

        pulse.data = struct;
        pulse.data.pulsetab = pulsetab;
        pulse.data.marktab = sortrows(marktab',1)';
        pulse.data.readout = readout;
        pulse.data.elem=pulsedef;  % Copy forward documentation.
        pulse.format = 'tab';
        
    otherwise
        error('Invalid format %s.', pulse.format);
end


        
