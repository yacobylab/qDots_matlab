function awglistgroups(name)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if nargin < 1
    pulsegroups = awgdata.pulsegroups;
else
    load([awgdata.grpdir, name], 'pulsegroups')
end

fprintf('Group\t%20s\tChannels\tLines\n', 'Name');
for l=1:length(pulsegroups)
    fprintf('%03d\t%25s\t%-5s\t%d\n', l, pulsegroups(l).name, sprintf('%d ',pulsegroups(l).chan(:)),...
            pulsegroups(l).npulse);
end
end
