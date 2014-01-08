function rinds = awglist(rng, name)
% inds = awglist(rng, name)
% rng: pulse indices to be shown. Default all if not given or empty.
% single number means last n pulses.
% name: show only pulses with this name.
% inds: indices of pulses displayed

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata

if ~isfield(awgdata, 'pulsedata')
    load(awgdata.datafile);
else
    pulsedata = awgdata.pulsedata;
end

if nargin < 1 || isempty(rng)
    rng = 1:length(pulsedata);
end

if rng(1) <= -1;
    rng = length(pulsedata)+rng:length(pulsedata);
end

pulses = pulsedata(rng);

if nargin >= 2
    inds = strmatch(name, {pulses.name});
    rng = rng(inds);
    pulses = pulses(inds);
end;

if nargout >= 1
    rinds = rng;
else
    fprintf('%-6s  %-10s  %-10s %-10s %-10s %-10s\n', 'pulse', 'name', 'clock','tau_rc','xdata','len')
    fprintf('--------------------------------------------------\n');
    for i = 1:length(rng);
        fprintf('%6d  %-10s  %5g  %10g %10g %f\n', rng(i), pulses(i).name, pulses(i).clk, pulses(i).taurc,pulses(i).xval,max(pulses(i).pulsetab(1,:)));
    end
end
