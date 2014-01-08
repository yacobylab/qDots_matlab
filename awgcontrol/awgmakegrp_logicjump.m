function awgmakegrp(pulsegroups, name)
% awgmakegrp(pulsegroups, name)
% pulsegroups is a struct array with fields
% pulses and nrep, ctrl (optional);
% ctrl: loop (does not work as it should)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if ~ isfield(pulsegroups, 'ctrl')
    [pulsegroups.ctrl] = deal('');
end

pls = [];
for i = 1:length(pulsegroups)
    pulsegroups(i).npulse = length(pulsegroups(i).pulses);
    pulsegroups(i).seqind = size(pls, 1)+1;
    if isempty(strfind(pulsegroups(i).ctrl, 'loop'))
        jumpline = -1;
    else
        jumpline = pulsegroups(i).seqind + 1;
    end
    if isfinite(pulsegroups(i).nrep)
        if length(pulsegroups(i).nrep) == 1
            pulsegroups(i).nrep = repmat(pulsegroups(i).nrep, 1, pulsegroups(i).npulse);
        end
        pulsegroups(i).npulse = length(pulsegroups(i).pulses);
        pls = [pls; [[awgdata.trigp; pulsegroups(i).pulses'; awgdata.offp], ...
            [1;  pulsegroups(i).nrep'; 0], zeros(pulsegroups(i).npulse+2, 2), ...
            [-ones(pulsegroups(i).npulse+1, 1); jumpline]]];
    else
        pls = [pls; [pulsegroups(i).pulses', zeros(pulsegroups(i).npulse, 3), ...
            [-ones(pulsegroups(i).npulse-1, 1); jumpline]]];
    end
end
awgmakeseq('', pls, [name, '.seq']);
save([awgdata.datadir, name], 'pulsegroups');

