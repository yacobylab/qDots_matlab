function pinf = awgmakeinf(pulsedef, pinf)
%pinf = awgmakeinf(pulsedef, pinf)
%
% pulsedef: struct array with fields type, time, val corresponding to pulse elements
%           to be concatenated.  type is a string specifying the type of element, 
%           val and time specicy pulse voltages and times, respectively.
%           Their meaning and the format of val depend on type.
% pinf: pulse description struct as passed to awgmakepulse. Pulsetab and marktab are set
%       according to pulsedef.
%
% Possible type strings and corresponding interpretation of val:
% raw: insert [time; val] into pulse table.
% fill: stretch this element to make the total pulse duration equal to time.
% wait: stay at val (row vector, one entry for each channel) for duration time.
% reload: relaod pulse at val (row vector, one entry for each channel).
%         time: [ramp time, wait time at load point, wait time at (0, 0) after load] 
% meas: measurement stage at [0, 0] for time(1), RF marker delayed by time(2) and
%       off time(3) before end of the stage.  [time(2) is lead delay,
%       time(3) is negative tail delay.
% ramp: ramp to val (row vector, one entry for each channel) in time.
% comp: measurement compensation at val(1:2) (one for each channel) for duration time(1). 
%       Ramps voltage to target and back over time(2) and time(3) at the beginning and 
%       end of the stage, respectively. If length(val)>=4, val(3:4) are used as final value.
%       The compensation value could be determined automatically, but this feature is not 
%       implemented yet.
% adprep: adiabatic ramp along second diagonal (epsilon) from val(1) to val(2), ramp duration time.
% adread: same, going the other way.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


pulsetab = zeros(3, 0);
marktab =  zeros(5, 0);
comppos = [];
fillpos = [];
readout = [];

for i = 1:length(pulsedef)
   
    switch pulsedef(i).type

        case 'raw'
            pulsetab = [pulsetab, [pulsedef(i).time; pulsedef(i).val]];
            
        case 'fill'
            fillpos = size(pulsetab, 2);
            filltime = pulsedef(i).time(1);
            
        case 'wait'
            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [1e-3, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
            pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(1:2)', 1, 2);
            
        case 'reload'
            pulsetab(1, end+(1:4)) = pulsetab(1, end) + cumsum(pulsedef(i).time([1 2 1 3]));
            pulsetab(2:3, end+(-3:0)) = [repmat(pulsedef(i).val(1:2)', 1, 2), zeros(2)];
            
        case 'meas'
            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [1e-3, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
            pulsetab(2:3, end+(-1:0)) = 0;
            marktab(:, end+1) = [pulsetab(1, end-2)+pulsedef(i).time(2); 0; 0; 0; pulsedef(i).time(1:3)*[1; -1; -1]];
            if ~isempty(pulsedef(i).val)
                readout(end+1, :) = [pulsedef(i).val, pulsetab(1, end-1), pulsedef(i).time(1)];
            end
                
        case 'ramp'
            pulsetab(1, end+1) = pulsetab(1, end) + pulsedef(i).time(1);
            pulsetab(2:3, end) = pulsedef(i).val(1:2);            
            
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
            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [1e-3, pulsedef(i).time(1)];
            pulsetab(2:3, end-1) = pulsedef(i).val(1)  * [-1 1];
            pulsetab(2:3, end) = pulsedef(i).val(2) * [-1 1];
            
        case 'adread'
            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [1e-3, pulsedef(i).time(1)];
            pulsetab(2:3, end-1) = pulsedef(i).val(2)  * [-1 1];
            pulsetab(2:3, end) = pulsedef(i).val(1)  * [-1 1];
            
        otherwise
            error('Invalid pulse element %i: %s.\n', i, pulsedef(i).type)           
    end
end

if ~isempty(comppos)
    pulsetab(2:3, comppos+(1:2)) = 2*repmat(compval(1:2)', 1, 2);
end

pulsetab(2:3, :) = pulsetab(2:3, :)./pinf.scale;
pinf = rmfield(pinf, 'scale');

if ~isempty(fillpos)
    filltime = filltime - pulsetab(1, end);
    if filltime < 0
        error('Pulse too long.');
    end
    pulsetab(1, fillpos+1:end) = pulsetab(1, fillpos+1:end) + filltime;
    
    mask = marktab(1, :) > pulsetab(1, fillpos);
    marktab(1, mask) = marktab(1, mask) + filltime; 
end

    

mask = all(abs(diff(pulsetab(2:3, :), [], 2)) < 1e-14);
pulsetab(:, [false, mask(2:end)&mask(1:end-1)]) = [];

pinf.pulsetab = pulsetab;
pinf.marktab = marktab;
