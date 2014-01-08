% Find the pulse line associated with a pulse group or pulse index.
% negative for groups, positive for pulse index.
% function seqind = awgseqind(pulses)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

function seqind = awgseqind(pulses,rep)
global awgdata;
if(isstruct(pulses))
    rep=[pulses.rep];
    pulses=[pulses.pulses];
end
seqind = nan(1, length(pulses));
for i = 1:length(pulses)
    if pulses(i) > 0
    if(exist('rep'))
        ind = find(pulses(i) == awgdata.seqpulses);
        ind=ind(rep(i));
    else        
        ind = find(pulses(i) == awgdata.seqpulses, 1);
    end
        if ~isempty(ind)
            seqind(i) = ind;
        end
    else
        seqind(i) = awgdata.pulsegroups(-pulses(i)).seqind;
    end
end
if any(isnan(seqind))
    fprintf('WARNING: Some pulses not present in sequence.\nHit Ctrl-C to abort, or any key to continue.\n');
    pause;
end
