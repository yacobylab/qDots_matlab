function awgsetamp(scale)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global tunedata;

inst = sminstlookup('AWG520');
v(1) = smquery(inst, 'SOUR1:VOLT:AMPL?', '%s', '%f');
v(2) = smquery(inst, 'SOUR2:VOLT:AMPL?', '%s', '%f');

v(3:4) = v(1:2) * 0.5 * scale./tunedata.plscal.fit(:, 6)';

logentry('Changed AWG amplitude from %.3f, %.3f V to %.3f, %.3f V.\n', v);

smprintf(inst, 'SOUR1:VOLT:AMPL? %f', v(3));
smprintf(inst, 'SOUR2:VOLT:AMPL? %f', v(4));

