
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


function polarize_B_variable(x, pulses, workp, duration, Bmeas)
% polarize for a time x(2) and wait for a time wait_time
% before measuring
%
% x: loop variables (fastest first)
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% duration: duration of polarization in sec.
% for this case put Brng into the second loop and take out the part with
% wait_time

smset({'GateR', 'GateL'}, workp(1:2));
smset('PulseLine', pulses(1)); %polarize
pause(duration);
smset('PulseLine', pulses(3)); %off
smset('B', Bmeas);
smset('PulseLine', pulses(2)); %measure
