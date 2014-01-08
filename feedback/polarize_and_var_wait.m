function polarize_and_variable_wait_fn(x, ctrl, pulses, workp, duration)
% polarize for a fixed time duration and wait for a time x(2)
% before measuring
%
% x: loop variables (fastest first)
% ctrl: not used
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



function variable_polarize_and_wait_fn(x, ctrl, pulses, workp, wait_time)
% polarize for a time x(2) and wait for a time wait_time
% before measuring
%
% x: loop variables (fastest first)
% ctrl: not used
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% wait_time: duration of waiting in seconds.

smset({'GateR', 'GateL'}, workp(1:2));
smset('PulseLine', pulses(1)); %polarize
pause(x(2));
smset('PulseLine', pulses(3)); %off
pause(wait_time);
smset('PulseLine', pulses(2)); %measure



function polarize_B_variable(x, ctrl, pulses, workp, duration)
% polarize for a time x(2) and wait for a time wait_time
% before measuring
%
% x: loop variables (fastest first)
% ctrl: not used
% pulses: 3-vector with pulse lines for polarization, measurement, and off.
% workp: GateR and GateL values to work at.
% duration: duration of polarization in sec.
% for this case put Brng into the second loop and take out the part with
% wait_time

smset({'GateR', 'GateL'}, workp(1:2));
smset('PulseLine', pulses(1)); %polarize
pause(duration);
smset('PulseLine', pulses(3)); %off
smset('B', 0.07);
smset('PulseLine', pulses(2)); %measure
