function awgcalcshift(plsnum)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if ~isfield(awgdata, 'pulsedata')
    load(awgdata.datafile);
else
    pulsedata = awgdata.pulsedata;
end

if nargin < 1
    plsnum = length(awgdata.shift)+1:length(awgdata.pulsedata);
end

for i = plsnum
    pulseinf = pulsedata(i);
    pulsetab = pulseinf.pulsetab;

    nchan = size(pulsetab, 1)-1;

    awgdata.shift(1:nchan, i) = 0.5 * sum((pulsetab(2:end, 2:end) + pulsetab(2:end, 1:end-1)).*...
        repmat((pulsetab(1, 2:end) - pulsetab(1, 1:end-1)), nchan, 1), 2)./(pulsetab(1, end) - pulsetab(1, 1));

end
