function polarizefnLR(x, ctrl, pulses, workp, duration, bfield)
% polarizefnLR(x, ctrl, pulses, workp, duration)
% Run polarization cycle with pulse pulses(1) at position workp for
% duration seconds.
%
% x: loop variables (fastest first)
% ctrl: [build | decay], [delay], [field].
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% duration: duration of polarization cycle in seconds.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;
    
if strfind(ctrl, 'build')
    dopol = x(2) > 60;
elseif strfind(ctrl, 'decay')       
    dopol = x(2) > 5 && x(2) <= 10; % measure base line, build up and decay
else        
    dopol = 1;
end

if strfind(ctrl, 'field');
    smset('B', bfield(1));
end

smset('PulseLine', pulses(2));
%pause(.1)
switch length(workp)
    case {3, 6}
        smset({'GateL', 'GateR', 'GateQPC'}, workp(1:3));
    case {2, 4}
        smset({'GateL', 'GateR'}, workp(1:2));
        
    case 9
        if dopol
            smset({'GateL', 'GateR', 'GateQPC'}, workp(1:3));
        else
            smset({'GateL', 'GateR', 'GateQPC'}, workp(7:9));
        end
        
end

if dopol  % measure base line, build up and decay
    smset('PulseLine', pulses(1));
else
    smset('PulseLine', pulses(3));
end

if strfind(ctrl, 'varpol')
    pause(x(3));
elseif strfind(ctrl, 'feedback')
    pause(fbdata.ctrlval(end));
else
    pause(duration);
end

if strfind(ctrl, 'delay')
    smset('PulseLine', pulses(3));
    pause(x(2)); % varying wait time after polarization
end  

if length(pulses) >= 2 && pulses(2) > 0
    smset('PulseLine', pulses(2));
end

switch length(workp)
    case 4
        smset({'GateL', 'GateR'}, workp(3:4));
    case {6 9}
        smset({'GateL', 'GateR', 'GateQPC'}, workp(4:6));
end

if strfind(ctrl, 'field');
    smset('B', bfield(2));
end

