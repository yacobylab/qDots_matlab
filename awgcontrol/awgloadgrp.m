function awgloadgrp(name, noload)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if isfield(awgdata, 'grpdir')
    load([awgdata.grpdir, name]);
else
    load([awgdata.datadir, name]);
end

awgdata.pulsegroups = pulsegroups;
infrep = false(1, length(pulsegroups));

for i = 1:length(pulsegroups)
    infrep(i) = ~isempty(pulsegroups(i).nrep) && ~all(isfinite(pulsegroups(i).nrep));
end


if isfield(awgdata, 'awg') && (nargin < 2 || noload == 0)       
    awgcntrl('stop');
    awgselect('seq', [name, '.seq']);    
    if isfield(pulsegroups, 'jump')
        awgcntrl('wait');
        err = query(awgdata.awg, 'SYST:ERR?');
        if ~strncmp(err, '0,"No error"', 12)
            error(err);
        end

        %start = [0 cumsum(2 * ~infrep + [pulsegroups.npulse])];
        % pulse and final off pulse added if all pulsereps are finite
        %start(end) = [];
        %start = start+~infrep; % add trigger pulses
        for i = 1:length(pulsegroups) % program jumps            
            for jump = pulsegroups(i).jump 
                fprintf(awgdata.awg, 'SEQ:ELEM%i:GOTO:IND %i', jump'+pulsegroups(i).seqind);
                fprintf(awgdata.awg, 'SEQ:ELEM%i:GOTO:STAT 1', jump(1)+pulsegroups(i).seqind);
            end
        end
        err = query(awgdata.awg, 'SYST:ERR?');
        if strncmp(err, '0,"No error"', 12)
            awgcntrl('start');            
        else
            error(err);
        end
    end
end

infrep = infrep(1:find(infrep, 1, 'last')); %avoid excess nans
npulse = sum([pulsegroups(1:length(infrep)).npulse])+2*sum(~infrep);

awgdata.seqpulses = nan(1, npulse);
for i = find(infrep);
    awgdata.seqpulses(pulsegroups(i).seqind +(0:pulsegroups(i).npulse-1)) = pulsegroups(i).pulses;
end
