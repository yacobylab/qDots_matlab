function polarizefn2LR(x, loop, pulsetab, workp)
% polarizefn2LR(x,  loop, pulsetab, workp, mode)
% 
% x: loop variables (fastest first)
% loop: index to loop variable to be used to determine pulse
% pulsetab: matrix with one polarization pulse per row
%           min loop val, max loop val, pulse, duration, workp ind, mode 
% min and max loop vals are inclucive. No pulse is set for pulse < 0.
% mode == 1: (default): run pulse for pulsetab(i, 4); 
% mode == 2: run pulse for x(pulsetab(i, 4));
% mode == 3: run pulse for pulsetab(i, 4) * max(fbdata.ctrlval(end)-0.5, .2)
% mode == 4: run pulse for pulsetab(i, 4) * max(0.5 - fbdata.ctrlval(end), .2)
% mode == 5: run pulse for fbdata.ctrlval(end)
% mode == 6: run pulse for pulsetab(i, 4) * abs(fbdata.ctrlval(end)) if
%            fbdata.ctrlval(end) > 0
% mode == 7: run pulse for pulsetab(i, 4) * abs(fbdata.ctrlval(end)) if
%            fbdata.ctrlval(end) < 0
% mode == 8: run pulse if listed in fbdata.pulseind

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


% Workp index <= 0 is ignored.


global fbdata;
pulses = find(pulsetab(:, 1) <= x(loop) & pulsetab(:, 2) >= x(loop))'; 

if isempty(pulses)
    pause(.1);
    return;
end

if size(pulsetab, 2) < 6
    pulsetab(:, 6) = 1;
end

gates = {'GateL', 'GateR', 'GateQPC'};
for i = pulses;
    
    if (pulsetab(i, 6) == 6 && fbdata.ctrlval(end) > 0) ...
            || (pulsetab(i, 6) == 7 && fbdata.ctrlval(end) < 0) ...
            || pulsetab(i, 6) == 8 && any(fbdata.pulseind == i)
        smset(gates(1:size(workp, 2)), workp(pulsetab(i, 5), :));
        smset('PulseLine', pulsetab(i, 3));
        if pulsetab(i, 4) > 0
            pause(pulsetab(i, 4) * abs(fbdata.ctrlval(end)));
        else
            pause(-pulsetab(i, 4));
        end
    elseif all(pulsetab(i, 6) ~= 6:8)
        if pulsetab(i, 5) > 0
            smset(gates(1:size(workp, 2)), workp(pulsetab(i, 5), :));
        end
        if pulsetab(i, 3) > 0
            smset('PulseLine', pulsetab(i, 3));
        end

        switch pulsetab(i, 6)
            case 1
                pause(pulsetab(i, 4));

            case 2
                pause(x(pulsetab(i, 4)));

            case 3
                pause(pulsetab(i, 4) * max(fbdata.ctrlval(end)-0.5, .2));

            case 4
                pause(pulsetab(i, 4) * max(0.5-fbdata.ctrlval(end), .2));

            case 5
                pause(fbdata.ctrlval(end));
        end
    end
    
end


