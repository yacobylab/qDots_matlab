function awgmakegrp(pulsegroups, name, chans)
% awgmakegrp(pulsegroups, name)
% pulsegroups is a struct array with fields
% pulses and nrep, ctrl (optional);
% ctrl: loop (only for first loop)
% nrep specifies how many times to repeat each pulse. It can be a single number 
% or a vector of the same length as pulses.
% If nrep is a single finite number, it is used for all but the last pulse, 
% which is run indefinitely (corresponding to nrep = 0).
% The exception are groups with only a single pulse, where nrep is not changed
% to 0.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if nargin < 3
    chans = awgdata.chans;
end

if ~ isfield(pulsegroups, 'ctrl')
    [pulsegroups.ctrl] = deal('');
end

pls = [];
for i = 1:length(pulsegroups)
    pulsegroups(i).npulse = length(pulsegroups(i).pulses);
    pulsegroups(i).seqind = size(pls, 1)+1;

    if isfinite(pulsegroups(i).nrep) % all finite
        if length(pulsegroups(i).nrep) == 1 && pulsegroups(i).npulse > 1
            %pulsegroups(i).nrep = repmat(pulsegroups(i).nrep, 1, pulsegroups(i).npulse);
            pulsegroups(i).nrep = [repmat(pulsegroups(i).nrep, 1, pulsegroups(i).npulse-1), 0];
        end
        pulsegroups(i).npulse = length(pulsegroups(i).pulses);
        if isempty(strfind(pulsegroups(i).ctrl, 'loop'))
            pls = [pls; [[awgdata.trigp; pulsegroups(i).pulses'; awgdata.offp], ...
                [1;  pulsegroups(i).nrep'; 0], zeros(pulsegroups(i).npulse+2, 2)]];
        else
            if i > 1
                error('Only first pulse group can be an endless loop.');
            end
            pls = [pls; [[pulsegroups(i).pulses'; awgdata.trigp], ...
                [pulsegroups(i).nrep'; 1], zeros(pulsegroups(i).npulse+1, 2)]];
            pls(end+(-1:0), end) = 1;
            pulsegroups(i).seqind = size(pls, 1);
        end
    elseif length(pulsegroups(i).nrep) == 1 
        pls = [pls; [pulsegroups(i).pulses', zeros(pulsegroups(i).npulse, 3)]];
    else
        pulsegroups(i).nrep(~isfinite(pulsegroups(i).nrep)) = 0;
        pls = [pls; [pulsegroups(i).pulses', pulsegroups(i).nrep', zeros(pulsegroups(i).npulse, 2)]];
    end
end
awgmakeseq('', pls, [name, '.seq'], chans);
if isfield(awgdata, 'grpdir')
    save([awgdata.grpdir, name], 'pulsegroups', 'chans');
else
    save([awgdata.datadir, name], 'pulsegroups', 'chans');
end


