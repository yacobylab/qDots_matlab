function pulseswitchfn(x)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;

persistent tstart;

if isempty(tstart)
    tstart = - inf(1, 2);
end

if strfind(fbdata.switchmode, 'rate')
    
    if length(fbdata.switchinds) < 2
        return;
    end

    if fbdata.fbval(end) < fbdata.switchlims(1) || (now - tstart(1)) * 24 * 3600 > fbdata.switchtime;
        fbdata.pulseind = fbdata.switchinds(:, 1)';
        tstart = [inf now];
    end

    if fbdata.fbval(end) > fbdata.switchlims(2) || (now - tstart(2)) * 24 * 3600 > fbdata.switchtime;
        fbdata.pulseind = fbdata.switchinds(:, 2)';
        tstart = [now inf];
    end

end

if ~isempty(strfind(fbdata.switchmode, 'check')) && fbdata.switchcount(1) == 0
    % check if still locked
    if fbdata.switchcount(2) > 0
        fbdata.switchcount(2) = fbdata.switchcount(2) - 1;
    elseif fbdata.fbval(end) < fbdata.switchlims(1)
        fbdata.switchcount = fbdata.switchcycles;
    elseif fbdata.fbval(end) > fbdata.switchlims(2)
        fbdata.switchcount = [-fbdata.switchcycles(1), fbdata.switchcycles(2)];
    end
end

if strfind(fbdata.switchmode, 'lock') % stay at different pulse for a number of cycles       
    if fbdata.switchcount(1) > 0
        fbdata.pulseind(2) = fbdata.switchinds(2);
        fbdata.switchcount(1) = fbdata.switchcount(1) - 1;
    elseif fbdata.switchcount(1) < 0
        fbdata.pulseind(2) = fbdata.switchinds(3);
        fbdata.switchcount(1) = fbdata.switchcount(1) + 1;
    else
        fbdata.pulseind(2) = fbdata.switchinds(1);
    end
end

if strfind(fbdata.switchmode, 'seq') % cycle through pulses?        
    ind = find(fbdata.pulseind == fbdata.switchinds);
    if  fbdata.fbval(end)*sign(fbdata.switchlims(ind)) > fbdata.switchlims(ind) || (now - tstart(1)) * 24 * 3600 > fbdata.switchtime(ind);
        tstart(1) = now;
        if ind == length(fbdata.switchinds)
            fbdata.pulseind = fbdata.switchinds(1);
        else
            fbdata.pulseind = fbdata.switchinds(ind+1);
        end
    end
end

if strfind(fbdata.switchmode, 'set') % move setpoint of HW pulses       
    if  fbdata.switchcount(1) > 0
        fbdata.switchcount(1) = fbdata.switchcount(1) - 1; % wait for Dbz to settle.
    elseif isfinite(fbdata.switchinds(1)) && fbdata.switchinds(1) ~= fbdata.pulseind;
        fbdata.pulseind = fbdata.pulseind + sign(fbdata.switchinds(1) - fbdata.pulseind);
        fbdata.switchcount(1) = fbdata.switchcount(2);
    end
end
