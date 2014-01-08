function scan = qpcfbconfig(scan)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


%scan.loops(2).procfn.dim = []; %scan.loops(1).npoints;
scan.loops(2).procfn.fn.fn = @qpcfeedbfn;
scan.loops(2).procfn.fn.args = {};
scan.loops(2).procfn(2).fn = []; % copy data

scan.loops(2).trafofn{end} = @qpctrafofn;
scan.loops(2).getchan{end+1} = scan.loops(2).setchan{end};

% fields of qpcdata:
% slp: [dVqpc/d(compensation gate),  dVqpc/d(scanned gate)]
% trco: [slope, offset] for transformation function.
% setp: target RF voltage.
% lim: lower and upper bound for Vqpc.

% reset log of gate voltages

function data = qpcfeedbfn(data)


global qpcdata;

m = mean(data);
qpcdata.trco(2) = qpcdata.trco(2) - (m - qpcdata.setp)/qpcdata.slp(1);

% keep track of last n compensation and gate values and readout values
% fit plane to update slopes if all valid, allow underrelaxation, set limits to slopes


function v = qpctrafofn(x, y)

global qpcdata;

v = max(min(x(2) * qpcdata.trco(1) + qpcdata.trco(2), qpcdata.lim(2)), qpcdata.lim(1)); 
