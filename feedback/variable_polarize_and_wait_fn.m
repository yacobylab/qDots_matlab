
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


function variable_polarize_and_wait_fn(x, pulses, workp, wait_time)
% polarize for a time x(2) and wait for a time wait_time
% before measuring
%
% x: loop variables (fastest first)
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% wait_time: duration of waiting in seconds.

smset({'GateR', 'GateL'}, workp(1:2));
smset('PulseLine', pulses(3)); %off
pause(wait_time);
smset('PulseLine', pulses(1)); %polarize
pause(x(2));
smset('PulseLine', pulses(2)); %measure


