function awgrmpulse(inds)
% delete pulse files and remove from database if at end.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

awglist(inds);

if ~isfield(awgdata, 'pulsedata')
    load(awgdata.datafile);
else
    pulsedata = awgdata.pulsedata;
end

if ~strcmp(input('Delete those pulses? (type yes)', 's'), 'yes')
    return
end

for i = sort(inds, 1, 'descend')
    name = pulsedata(i).name;

    if ~isempty(name)
        name = ['_', name];
    end

    if isfield(awgdata, 'plsgen');
        name = [awgdata.plsgen, name];
    end

    for j = 1:size(pulsedata(i).pulsetab, 1)-1;        
        delete(sprintf('%sp%06d%s_%1d.wfm', awgdata.datadir, i, name, j));
    end
    
    if i==length(pulsedata)
        pulsedata(i) = [];
    end
end


if ~isfield(awgdata, 'pulsedata')
    save(awgdata.datafile, 'pulsedata');
else
    awgdata.pulsedata = pulsedata;
end
