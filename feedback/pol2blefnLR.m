function pol2blefnLR(x, ctrl, pulses, workp, duration, nswitch)
% pol2blefnlR(x, ctrl, pulses, workp, duration, nswitch)
% Simplified polarizefn designed to switch pulses and duration 
% part way through a measurement.
% duration seconds.
%
% x: loop variables (fastest first)
% ctrl: unused
% pulses: 3-vector with pulse lines for twp polarizations and off.
% workp: GateR and GateL values to work at, 3x3 matrix
% duration: duration of polarization cycle in seconds, one for each pulse

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.




smset('PulseLine', pulses(3));

ind = 2-(x(3)<=nswitch);

smset({'GateL', 'GateR', 'GateQPC'}, workp(ind, :));
smset('PulseLine', pulses(ind));
pause(duration(ind));
smset('PulseLine', pulses(3));
smset({'GateL', 'GateR', 'GateQPC'}, workp(3, :));


