function polarizefn2bLR(x, loop, pulsetab, workp)
% polarizefn2bLR(x,  loop, pulsetab, workp, mode)
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

%workp(:, 1:length(fbdata.offset)) = workp(:, 1:length(fbdata.offset))...
%    + repmat(min(max(fbdata.offset, -5e-3),  5e-3), 1, size(workp, 1));
%workp(:, 2) = workp(:, 2)+ min(max(fbdata.offset, -5e-3),  5e-3);
if size(fbdata.offset, 1) == size(workp, 1)
    workp(:, 1:2) = workp(:, 1:2) + min(max(fbdata.offset, -5e-3),  5e-3);
else
    workp(:, 1:2) = workp(:, 1:2) + repmat(min(max(fbdata.offset, -5e-3),  5e-3), size(workp, 1), 1);
end
sweep = min(max(fbdata.sweep, -5e-3), 5e-3);


gates = {'GateL', 'GateR', 'GateQPC'};
for i = pulses;
       
    switch pulsetab(i, 6)
        case 0
            dt = 0;
        case 1
            dt = pulsetab(i, 4);
            
        case 2
            dt = x(pulsetab(i, 4));
            
        case 3
            dt =pulsetab(i, 4) * max(fbdata.ctrlval(end)-0.5, .2);
            
        case 4
            dt = pulsetab(i, 4) * max(0.5-fbdata.ctrlval(end), .2);
            
        case 5
            dt = fbdata.ctrlval(end);
            
        case 6
            if fbdata.ctrlval(end) <= 0; continue; end
            
        case 7
            if fbdata.ctrlval(end) >= 0; continue; end

        case 8            
            if all(fbdata.pulseind ~= i); continue; end
            
    end
    
    if any(pulsetab(i, 6) == 6:8)
        if pulsetab(i, 4) > 0
            dt = pulsetab(i, 4) * abs(fbdata.ctrlval(end));
        else
            dt = -pulsetab(i, 4);
        end
    end   

    if pulsetab(i, 5) > 0
        smset(gates(1:size(workp, 2)), workp(pulsetab(i, 5), :));
    end
    
    if pulsetab(i, 3) <= 0 % do nothing if pulse invalid
        continue; % does not allow pauses without pulse, but that should be fine
    end

    if abs(sweep) > 0 && pulsetab(i, 5) > 0
        smset(gates(2), workp(pulsetab(i, 5), 2) + sweep, -dt/abs(sweep));
        smset('PulseLine', pulsetab(i, 3));
        smatrigfn(smchaninst(gates(2)));  %obsolete if HW triggering      
    else
        smset('PulseLine', pulsetab(i, 3));
    end
    pause(dt);
end



