function scan = smasweepHP8350(scan)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


inst = sminstlookup('HP8350B');

f = smquery(inst, 'OP CW', '%s\n', '%f');

smprintf(inst, sprintf('FA %f HZ', scan.loops(1).rng(1)));
smprintf(inst, sprintf('FB %f HZ', scan.loops(1).rng(2)));
smprintf(inst, sprintf('ST %f SC', scan.loops(1).npoints*abs(scan.loops(1).ramptime)*.9));
smprintf(inst, 'T3');

% check for other cleanupfn?
scan.cleanupfn.fn = @smaconfigwrap;
scan.cleanupfn.args = {@smprintf, inst, sprintf('CW %f HZ', f)};

% set trigfn?
