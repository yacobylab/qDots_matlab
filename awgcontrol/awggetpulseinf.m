function pulse = awggetpulseinf(inds)
% pulse = awggetpulseinf(inds)
% returnd pulse definitions for pulses whose indices are ginven in inds.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if isfield(awgdata, 'pulsedata')
    pulse = awgdata.pulsedata(inds);
else
    load(awgdata.datafile);
    pulse = pulsedata(inds);
end
