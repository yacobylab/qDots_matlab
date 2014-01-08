function awgsavedata
% awgsavedata
% save awgdata in plsdata.grpdir, with name generated from date and time.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global awgdata;
global plsdata;

data = rmfield(awgdata, 'awg');
time = clock;
save(sprintf('%sawgdata_%02d%02d%02d_%02d%02d', plsdata.grpdir, mod(time(1), 100), time(2:5)), 'data','-v6');
