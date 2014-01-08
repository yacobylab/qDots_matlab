function awgadd(groups)
% awgadd(groups)
% Add groups to end of sequence. Store group name and target index in
% awgdata.pulsgroups.name, seqind.
% For sequence combined groups (only), pulseind is a 2D array.  Each
% row corresponds to a groups, each column to a varpar in this group.
% ie, pulseind = [ 1 1 1 1 ; 1 2 3 4] will use pls 1 from group 1, pulse
% 1-4 on group 2.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global plsdata;
global awgdata;
astart=toc;
awgcntrl('clr');
awgcntrl('stop');

if ~iscell(groups)
    groups = {groups};
end

dosave = false; % keeps track of whether awgdata changed.
gstart=toc;
for k = 1:length(groups)
    load([plsdata.grpdir, 'pg_', groups{k}]);
    
    while plsinfo('stale',groups{k})
        if ~awgdata(1).quiet
          fprintf('Latest pulses of group %s not loaded; %s > %s.\n', groups{k}, ...
            datestr(lastupdate), datestr(plslog(end).time(end)));
        end
        tstart=toc;
        plsmakegrp(groups{k},'upload');
        ts2 = toc;
        awgcntrl('wait');
        %  fprintf('Load time=%f secs, wait time=%f\n',toc-tstart,toc-ts2);
        load([plsdata.grpdir, 'pg_', groups{k} '.mat']);
    end
    
    
    if strcmp(grpdef.ctrl(1:min([end find(grpdef.ctrl == ' ', 1)-1])), 'grp')...
            && ~isempty(strfind(grpdef.ctrl, 'seq')) % combine groups at sequence level.
        
        % retrieve channels of component groups
        clear chan;
        lastload = 0;
        for m = 1:length(grpdef.pulses.groups)
            gd=plsinfo('gd', grpdef.pulses.groups{m});
            pl=plsinfo('log', grpdef.pulses.groups{m});
            if pl(end).time > lastload
                lastload = pl(end).time;
            end
            rf={'varpar'}; % Required fields that may be missing
            for qq=1:length(rf) % Put in empty copies of all required fields.
                if ~isfield(gd,rf{qq})
                    gd=setfield(gd,rf{qq},[]);
                end
            end
            chan(m) = orderfields(gd);
        end
        chan = {chan.chan};
        seqmerge = true;
    else
        if ~isfield(grpdef, 'pulseind')
            grpdef.pulseind = 1:size(zerolen, 1);
        end
        
        seqmerge = false;
    end
    
    if ~isfield(grpdef, 'nrep')
        grpdef.nrep = 1;
    end
    
    for a=1:length(awgdata)
        npls = size(grpdef.pulseind, 2);
        nchan = length(awgdata(a).chans); % alternatively use awgdata or data size
        usetrig = (grpdef.nrep(1) ~= Inf) && isempty(strfind(grpdef.ctrl, 'notrig'));
        
        if isempty(awgdata(a).pulsegroups)
            startline = 1;
            gind = [];
        else
            gind = strmatch(grpdef.name, {awgdata(a).pulsegroups.name}, 'exact');
            if ~isempty(gind)
                startline = awgdata(a).pulsegroups(gind(1)).seqind;
                if npls + usetrig ~= sum(awgdata(a).pulsegroups(gind(1)).npulse);
                    error('Number of pulses changed in group %s. Use awgrm first!', grpdef.name);
                end
            else
                startline = awgdata(a).pulsegroups(end).seqind + sum(awgdata(a).pulsegroups(end).nline);
            end
        end
        
        if ~exist('zerolen','var')  % For sequence combined groups, getting this is hard.
              zerolen = plsinfo('zl',grpdef.name);
        end
        if isempty(gind) % group not loaded yet, extend sequence
            
            gind = length(awgdata(a).pulsegroups)+1;
            awgdata(a).pulsegroups(gind).name = grpdef.name;
            awgdata(a).pulsegroups(gind).seqind = startline;
            awgdata(a).pulsegroups(gind).lastupdate = lastupdate;
            awgdata(a).pulsegroups(gind).npulse = [npls usetrig];
            if strfind(grpdef.ctrl,'pack')
                awgdata(a).pulsegroups(gind).nline = 1+usetrig;
                % Hack alert; way too much code assumes zl == pulselen.  For
                % packed groups, we ignore it and work out the correct length
                % ourselves.
                zlmult=npls;
                npls=1;
            else
                zlmult=1;
                awgdata(a).pulsegroups(gind).nline = npls+usetrig;
            end
            fprintf(awgdata(a).awg, sprintf('SEQ:LENG %d', startline + awgdata(a).pulsegroups(gind).nline-1));
            dosave = 1;
        else
            if strfind(grpdef.ctrl,'pack')
                zlmult = npls;
                npls=1;
            else
                zlmult=1;
            end
            if any(awgdata(a).pulsegroups(gind).nrep ~= grpdef.nrep) || ...
                    (isfield(plslog(end),'readout') && isfield(plslog(end).readout,'readout') && (any(any(awgdata(a).pulsegroups(gind).readout.readout ~= plslog(end).readout.readout)))) || ...
                    any(any(awgdata(a).pulsegroups(gind).zerolen ~= zerolen)) % nrep or similar changed
                dosave = 1;
            end
        end
        
        if ~isfield(grpdef, 'jump')
            if strfind(grpdef.ctrl, 'loop')
                grpdef.jump = [npls; 1];
            else
                grpdef.jump = [];
            end
        end
        
        
        % Cache frequently used information in awgdata
        awgdata(a).pulsegroups(gind).nrep = grpdef.nrep;
        if seqmerge
            awgdata(a).pulsegroups(gind).lastload = lastload;
        else
            awgdata(a).pulsegroups(gind).lastload = plslog(end).time(1);
        end
    
        awgdata(a).pulsegroups(gind).zerolen = zerolen; 
        if isfield(plslog(end),'readout')
          awgdata(a).pulsegroups(gind).readout=plslog(end).readout;
        else
          awgdata(a).pulsegroups(gind).readout = plsinfo('ro',grpdef.name);
        end
        
        if usetrig            
            for j = 1:nchan
                if mod(j,2) == 1
                  fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "trig_%08d"', startline, j, awgdata(a).triglen));
                else
                  fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', startline, j, awgdata(a).triglen, awgdata(a).zerochan(j)));
                end
            end
            if isfield(awgdata(a),'slave') && ~isempty(awgdata(a).slave) && (awgdata(a).slave)
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:TWAIT 1\n', startline));
            end
        end
        
        
        for i = 1:npls
            ind = i-1 + startline + usetrig;
            if ~seqmerge % pulses combined here.
                for j = 1:nchan
                    ch = find(awgdata(a).chans(j) == grpdef.chan);
                    if ~isempty(ch) &&  zerolen(grpdef.pulseind(i), ch) < 0
                        % channel in group and not zero
                        fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, j, ...
                            grpdef.name, grpdef.pulseind(i), ch));
                    else
                        % hack alert. We should really make zerolen a cell array.  fixme.
                        fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', ind, j, ...
                            zlmult*abs(zerolen(grpdef.pulseind(i), 1))*awgdata(a).clk/awgdata(1).clk,awgdata(a).zerochan(j)));
                    end
                end
            else
                for m = 1:length(grpdef.pulses.groups)
                    for j = 1:length(chan{m}) % channels of component groups
                        ch = find(awgdata(a).chans == chan{m}(j));
                        if ~isempty(ch)
                            %if 1 % zero replacement not implemented
                            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, ch, ...
                                grpdef.pulses.groups{m}, grpdef.pulseind(m, i), ch));
                        else
                            error('This won''t work');
                        end
                        %else
                        %fprintf(awgdata.awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d"', ind, awgdata.chans(j), ...
                        %    abs(zerolen(grpdef.pulseind(i), 1))));
                        %end
                    end
                end
            end
            if ~isempty(strfind(grpdef.ctrl,'single'))
                if ~isempty(strfind(grpdef.ctrl,'loop'))
                   error('trying to use ''loop'' and ''single'' for same group') 
                end
                %for a single, we want the group to run just once, then go
                %to pulseline 1. this will set it to run once
                fprintf(awgdata(a).awg, 'SEQ:ELEM%d:LOOP:INF 0', ind); % default
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:COUN 1', ind));
            elseif grpdef.nrep(min(i, end)) == Inf  || grpdef.nrep(min(i, end)) == 0 ...
                    || (i == npls && isempty(strfind(grpdef.ctrl, 'loop')) && (isempty(grpdef.jump) || all(grpdef.jump(1, :) ~= i)))
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:INF 1', ind));
            else
                fprintf(awgdata(a).awg, 'SEQ:ELEM%d:LOOP:INF 0', ind); % default
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:COUN %d', ind, grpdef.nrep(min(i, end))));
            end
            
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 0', ind));
            
            if grpdef.nrep(min(i, end)) == Inf && isreal(grpdef.pulses) &&  ...
                    (length(awgdata(a).seqpulses) < ind || awgdata(a).seqpulses(ind) ~= grpdef.pulses(grpdef.pulseind(i)));
                dosave = 1;
                awgdata(a).seqpulses(ind) = grpdef.pulses(grpdef.pulseind(i));
            end
            if ~mod(i, 100) && ~awgdata(1).quiet
                fprintf('%i/%i pulses added.\n', i, npls);
            end
        end
        %fprintf('Group load time: %g secs\n',toc-gstart);
        
        jstart=toc;
        % event jumps
        %SEQ:ELEM%d:JTARget:IND
        %SEQ:ELEM%d:JTARget:TYPE
        
        if ~isempty(strfind(grpdef.ctrl,'single'))
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:IND 1', ind)); %hack, that ind is the end of the sequence, left over from the prvious for loop!
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 1', ind));
        end
        
        for j = 1:size(grpdef.jump, 2)
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:IND %d', startline+usetrig-1 + grpdef.jump(:, j)));
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 1', startline+usetrig-1 + grpdef.jump(1, j)));
        end
    end
    if ~exist('seqlog','var')
        seqlog.time = now;
    else
        seqlog(end+1).time = now;
    end
    seqlog(end).nrep = grpdef.nrep;
    seqlog(end).jump = grpdef.jump;
    
    save([plsdata.grpdir, 'pg_', groups{k}], '-append', 'seqlog');
    %fprintf('Jump program time: %f secs\n',toc-jstart);
    wstart=toc;
    awgcntrl('wait');
    %fprintf('Wait time: %f secs; total time %f secs\n',toc-wstart,toc-astart);
    nerr=0;
    for a=1:length(awgdata(a))
        err=query(awgdata(a).awg, 'SYST:ERR?');
        if ~isempty(strfind(err, 'No error'))
            nerr=nerr+1;
        end
    end
    if nerr == 0 
        if ~awgdata(1).quiet
          fprintf('Added group %s on index %i. %s', grpdef.name, gind, err);
        end
        logentry('Added group %s on index %i.', grpdef.name, gind);
    end
end
if dosave
    awgsavedata;
end
