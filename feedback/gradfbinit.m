function gradfbinit

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;

fbdata.setp = 0;
fbdata.fbon = 0;
fbdata.plot = 1;
fbdata.igain  = .005;
fbdata.pgain  = .2;
fbdata.method = 3;

npoints = 200;
fbdata.time = nan(1, npoints);
fbdata.fbval = nan(1, npoints);
fbdata.intval = nan(1, npoints);
fbdata.ctrlval = nan(1, npoints);
fbdata.intval(end) = 0;
fbdata.ctrlval(end) = .1;

fbdata.dataind = 1;
fbdata.pulses = 1:32;
fbdata.freq = 1/10;
