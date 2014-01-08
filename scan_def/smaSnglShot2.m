function scan = smaSnglShot2(scan, procpar, ctrl, samplerate, nrep)
% smaSnglShot(scan, procpar, ctrl, samplerate, nrep)
% ctrl: dec: driver level averaging and decimation
%       raw: save raw data.
% samplerate is optional, actual rate or actual and hw rate.
% defaults are taken from configfn and consts.
% nrep: number of outer repetitions (optional)
%       Note: correlation functions are ill defined if nrep > 1 and each
%           pulse is repeated.
% if 'dec' is given, samplerate(1) is ignored.
% samplerate(2) is  taken from consts or configfn if not given.
%
% Assumes scan configured with smabufconfig2 and fast mode set in confSeq.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;
global plsdata;
global smdata;

if nargin < 3
    ctrl = '';
end

if nargin < 5 
    nrep = 1;
end

%find loop?
loop = 1;
scan.loops(loop).procfn = struct;


ic = smchaninst(scan.loops(loop).getchan);
nchan = sum(ic(:, 1) == ic(1, 1)); % # channels from same instrument as first

plsgrp = scan.data.pulsegroups(1);

if isfield(scan, 'consts')
    srch = strmatch('samprate', {scan.consts.setchan});
    if isempty(srch)
        srch = length(scan.consts) + 1;
    end
else
    srch = 1;
end

 
procpar.readout = [];
if strfind(ctrl, 'dec')   % using driver level downsampling
    if nargin < 4 || length(samplerate) < 2
        if ~isfield(scan, 'consts') || srch > length(scan.consts)
            error('Samplerate not given.');
        end           
        samplerate(2) = scan.consts(srch).val;
    end

    %readout = plsinfo('ro', plsgrp.name);
    % see smaSnglShot.m for mask settign code - obsoloete, mask set by configfn.
    readout = 1;  % FIXME; the meaning of readouts has changed a bit :(
    procpar.n = [size(readout, 1), plsgrp.nrep(1), plsgrp.npulse(1)];
    nmeas = nchan * procpar.n(3) * procpar.n(1);
else % configure readout sections
    if nargin < 4 || isempty(samplerate)
        samplerate = scan.configfn(1).args{3}(2);
    else
        scan.configfn(1).args{3}(2) = samplerate(1);
    end

    tpls = 0; % allow offset!
    %     for i = 1:plsgrp.npulse % add readouts for all pulses.
    %         readout = vertcat(awgdata.pulsedata(plsgrp.pulses(i)).readout);
    %         if isempty(readout)
    %             tpls = tpls + awgdata.pulsedata(plsgrp.pulses(i)).pulsetab(1, end);
    %         else
    %             readout(:, 2) = readout(:, 2) + tpls;
    %             tpls = tpls + awgdata.pulsedata(plsgrp.pulses(i)).pulsetab(1, end);
    %             procpar.readout = [procpar.readout, round(readout(:, 2:3)'*samplerate(1)*pls1.tbase/pls1.clk)];
    %         end
    %     end

    if isfield(procpar, 'avtime');
        ntime = length(procpar.avtime);
    else
        ntime = 0;
    end
    
    zl = plsinfo('zl', plsgrp.name);
    readout2 = plsinfo('ro', plsgrp.name);

    for i = 1:plsgrp.npulse(1) % add readouts for all pulses.
        %readout = repmat(awgdata.pulsedata(plsgrp.pulses(i)).readout, max(1, ntime), 1);
        readout = repmat(readout2, max(1, ntime), 1);
        % should use pulse specific readout here. Channel specific readout not supported so far
        if ~isempty(readout)
            if ntime > 0
                nread = size(awgdata.pulsedata(plsgrp.pulses(i)).readout, 1);
                readout(:, 3) = reshape(repmat(procpar.avtime, nread, 1), nread * ntime, 1);
            end
            readout(:, 2) = readout(:, 2) + tpls;
            procpar.readout = [procpar.readout, round(readout(:, 2:3)'*samplerate(1)*plsdata.tbase/awgdata.clk)];
        end

        tpls = tpls + abs(zl(min(i, end), 1))./plsdata.tbase; 
        % taking min is a dirty bug fix for composite groups.
    end

    % would need only single pulse if expanding
    % procpar.readout = round([awgdata.pulsedata(plsgrp.pulses).readout]*samplerate)
    procpar.n = [round(abs(zl(1, 1))*samplerate(1)/awgdata.clk), plsgrp.nrep(1), plsgrp.npulse(1)];   

    nmeas = nchan*size(procpar.readout, 2);
    %scan.configfn(2).args{3} = []; % removed 01/28. Not clear what it was for.
end

if nrep > 1
    procpar.n(4) = nrep;
end

if length(samplerate) >= 2
    scan.consts(srch).setchan = 'samprate';
    scan.consts(srch).val = samplerate(2);
    smset('samprate', samplerate(2)); % ensure check below is correct
end

scan.configfn(1).args{3}(1) = prod(procpar.n);


%run configfn for testing
% smabufconfig does not return values. Need to change that to perform this check.
% [n(1), n(2)] = scan.configfn.args{1}(scan.configfn.args{2:end});
% if any(n ~= [scan.configfn.args{3:4}])
%     error('Invalid sampling parameter.');
% end

% assumes all pulsegroups have the same pattern and all pulses in a group are repeated the same number
% of times and have the same length. Readout positions could in principle vary.
% May also need to allow a sequence of different single pulses. In that case, npulse = 1, 
% nsamples = Sum over all pulses, nrep determined from record length (not specified in group).
% code for readout above should work the same way.
 

if ischar(scan.loops(loop).getchan)
    inchan = 1;
    npc = 0;
else
    npc = length(scan.loops(loop).getchan)-nchan;
    inchan = 1:nchan;
end
if strfind(ctrl, 'raw')    
    npc = npc+nchan;
    coff=0;
else
    coff=nchan;
end
% npc = #unprocessed channels


ndata = length(procpar.datadef);
[scan.loops(loop).procfn(1:npc+ndata).fn] = deal([]);

% copy aux channels
if ndata > nchan % first copy non-processed channels back
    scan.loops(loop).procfn(1).fn.fn = [];
    scan.loops(loop).procfn(1).fn.inchan = coff + (1:npc);
    scan.loops(loop).procfn(1).fn.outchan = ndata + (1:npc);
    
    %set datadim for copied channels
    for i = 1:npc
        %inst =
        %smdata.channels(smchanlookup(scan.loops(loop).getchan(i+nchan))).instchan;
        inst = smdata.channels(smchanlookup(scan.loops(loop).getchan(i+coff))).instchan;
        if isempty(smdata.inst(inst(1)).datadim)
            scan.loops(loop).procfn(i+ndata).dim = 1;
        else
            scan.loops(loop).procfn(i+ndata).dim = smdata.inst(inst(1)).datadim(inst(2), 1); %only vectors implemented. OK with 0?
        end
    end
end


% alternative putting additional channels at end - needs modification 
%  of smrun, see software.txt, 12/28/09
%     scan.loops(loop).procfn(1).fn = pf;
%     [scan.loops(loop).procfn(2:ndata).fn] = deal([]);
% 
%     for i = 1:length(scan.loops(loop).getchan)-ndata
%         scan.loops(loop).procfn(ndata+i).fn.fn = '=';
%         scan.loops(loop).procfn(data+i).fn.inchan = nchan+i;
%     end

% set data dimensions for each channel
for i = 1:ndata    
    switch procpar.datadef(i).type

        case {'ave', 'gave'}
            if strcmp(procpar.datadef(i).par, 'readout')
                procpar.datadef(i).par = {procpar.readout(1, 1)+(1:procpar.readout(2, 1))};
            end
            if isempty(procpar.datadef(i).par)
                scan.loops(loop).procfn(i).dim = [nchan, procpar.n(1)];
            else
                scan.loops(loop).procfn(i).dim = [nchan, length(procpar.datadef(i).par{1})];
            end
            
            if strcmp(procpar.datadef(i).type, 'ave')
                scan.loops(loop).procfn(i).dim(1) = procpar.n(3) * scan.loops(loop).procfn(i).dim(1);
            end
            
        %case 'gave'
        %    scan.loops(loop).procfn(i).dim = [nchan, procpar.n(1)];

        case {'mean', 'bin'}
            if length(procpar.datadef(i).par) >= 1
                n = length(procpar.datadef(i).par{2});
            else
                n = nmeas;
            end
            scan.loops(loop).procfn(i).dim = [n, procpar.n(2)*nrep];

        case 'hist'
            if length(procpar.datadef(i).par) >= 2
                n = length(procpar.datadef(i).par{2});
            else
                n = nmeas;
            end
            scan.loops(loop).procfn(i).dim = [n, size(procpar.datadef(i).par{1}, 2)-1];

        case 'ghist'
            scan.loops(loop).procfn(i).dim = size(procpar.datadef(i).par{1}, 2)-1;

        case {'mom', 'bmom'}
            if isempty(procpar.datadef(i).par{1})
                scan.loops(loop).procfn(i).dim = nmeas;
            else
                scan.loops(loop).procfn(i).dim = size(procpar.datadef(i).par{1}, 2);
            end

        case {'corr', 'bcorr'}
            if isempty(procpar.datadef(i).par{1})
                n = nmeas;
            else
                n = size(procpar.datadef(i).par{1}, 2);
            end
            scan.loops(loop).procfn(i).dim = [n, 2*procpar.datadef(i).par{2}+1];
            
    end
end

% all processing happens on new data, no cumulative summing implemented.
scan.loops(loop).procfn(1).fn(2).fn = @smpSnglShot;
scan.loops(loop).procfn(1).fn(2).outchan = 1:ndata;
scan.loops(loop).procfn(1).fn(2).args = {procpar};
scan.loops(loop).procfn(1).fn(2).inchan = inchan;

