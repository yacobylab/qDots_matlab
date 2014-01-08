function polarize_and_variable_wait_fn(x, pulses, workp, duration)
% polarize for a fixed time duration and wait for a time x(2)
% before measuring
%
% x: loop variables (fastest first)
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% duration: duration of polarization cycle in seconds.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.



smset({'GateR', 'GateL'}, workp(1:2));
smset('PulseLine', pulses(1)); %polarize
pause(duration);
smset('PulseLine', pulses(3)); %off
pause(x(2));
smset('PulseLine', pulses(2)); %measure

