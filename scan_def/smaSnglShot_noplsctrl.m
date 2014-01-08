function scan = smaSnglShot(scan, procpar, ctrl, samplerate, nrep)
% smaSnglShot(scan, procpar, ctrl, samplerate, nrep)
% ctrl: dec: driver level averaging and decimation
%       raw: save raw data.
% samplerate is optional, actual rate or actual and hw rate.
% defaults are taken from configfn and consts.
% nrep: number of outer repetitions (optional)
%       Note: correlation functions are ill defined if nrep > 1 and each
%           pulse is repeated.
% if 'dec' is given, samplerate(1) is ignored.
% samplerate(2) is  taken from consts if not given.
%
% change scan.data.plsgrp.npulse to record not all pulses.
global awgdata;
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
nchan = sum(ic(:, 1) == ic(end, 1)); % # channels from same instrument as last
%nchan = length(scan.loops(loop).getchan); % allow other channels?

plsgrp = scan.data.pulsegroups(1);
pls1 = awgdata.pulsedata(plsgrp.pulses(1));


if isfield(scan, 'consts')
    srch = strmatch('samprate', {scan.consts.setchan});
    if isempty(srch)
        srch = length(scan.consts) + 1;
    end
else
    srch = 1;
end

 
procpar.readout = [];
if strfind(ctrl, 'dec')   %
    if nargin < 4 || length(samplerate) < 2
        if ~isfield(scan, 'consts') || srch > length(scan.consts)
            error('Samplerate not given.');
        end           
        samplerate(2) = scan.consts(srch).val;
    end

    readout = round(pls1.readout(1, 2:3) * samplerate(2)*pls1.tbase/pls1.clk);
    n = pls1.pulsetab(1, end) * samplerate(2)*pls1.tbase/pls1.clk;
    smdata.inst(ic(end, 1)).data.mask = false(1, n);
    smdata.inst(ic(end, 1)).data.mask(readout(1)+1:readout(1)+readout(2)) = true;
    procpar.n = [size(pls1.readout, 1), plsgrp.nrep(1), plsgrp.npulse];
    
    nmeas = nchan * procpar.n(3) * procpar.n(1);
    scan.configfn.args{4} = samplerate(2)/n;
else
    if nargin < 4
        samplerate = scan.configfn.args{4};
    else
        scan.configfn.args{4} = samplerate(1);
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
    
    for i = 1:plsgrp.npulse % add readouts for all pulses.
        readout = repmat(awgdata.pulsedata(plsgrp.pulses(i)).readout, max(1, ntime), 1);
    
        if ~isempty(readout)
            if ntime > 0
                nread = size(awgdata.pulsedata(plsgrp.pulses(i)).readout, 1);
                readout(:, 3) = reshape(repmat(procpar.avtime, nread, 1), nread * ntime, 1);
            end
            readout(:, 2) = readout(:, 2) + tpls;
            procpar.readout = [procpar.readout, round(readout(:, 2:3)'*samplerate(1)*pls1.tbase/pls1.clk)];
        end

        tpls = tpls + awgdata.pulsedata(plsgrp.pulses(i)).pulsetab(1, end);
    end


    % would need only single pulse if expanding
    % procpar.readout = round([awgdata.pulsedata(plsgrp.pulses).readout]*samplerate)
    procpar.n = [round(pls1.pulsetab(1, end)*samplerate(1)*pls1.tbase/pls1.clk), ...
        plsgrp.nrep(1), plsgrp.npulse];   

    nmeas = nchan*size(procpar.readout, 2);
    smdata.inst(ic(end, 1)).data.mask = [];
end

if nrep > 1
    procpar.n(4) = nrep;
end

if length(samplerate) >= 2
    scan.consts(srch).setchan = 'samprate';
    scan.consts(srch).val = samplerate(2);
    smset('samprate', samplerate(2)); % ensure check below is correct
end

scan.configfn.args{3} = prod(procpar.n);


%run configfn for testing
[n(1), n(2)] = scan.configfn.args{1}(scan.configfn.args{2:end});
if any(n ~= [scan.configfn.args{3:4}])
    error('Invalid sampling parameter.');
end

% assumes all pulsegroups have the same pattern and all pulses in a group are repeated the same number
% of times and have the same length. Readout positions could in principle vary.
% May also need to allow a sequence of different single pulses. In that case, npulse = 1, 
% nsamples = Sum over all pulses, nrep determined from record length (not specified in group).
% code for readout above should work the same way.
 

if ischar(scan.loops(loop).getchan)
    pf.inchan = 1;
    npc = 0;
else
    npc = length(scan.loops(loop).getchan)-nchan;
    pf.inchan = npc+1:npc+nchan; 
end
if strfind(ctrl, 'raw')    
    npc = npc+nchan;
end
pf.fn = @smpSnglShot;
ndata = length(procpar.datadef);
pf.outchan = npc+1:npc+ndata;

pf.args = {procpar};

[scan.loops(loop).procfn(1:npc+ndata).fn] = deal([]);
scan.loops(loop).procfn(npc+1).fn = pf;

% alternative putting additional channels at end - needs modification 
%  of smrun, see software.txt, 12/28/09
%     scan.loops(loop).procfn(1).fn = pf;
%     [scan.loops(loop).procfn(2:ndata).fn] = deal([]);
% 
%     for i = 1:length(scan.loops(loop).getchan)-ndata
%         scan.loops(loop).procfn(ndata+i).fn.fn = '=';
%         scan.loops(loop).procfn(data+i).fn.inchan = nchan+i;
%     end

for i = 1:ndata    
    switch procpar.datadef(i).type

        case 'ave'
            scan.loops(loop).procfn(i+npc).dim = procpar.n(1) * procpar.n(3)* nchan;
            
        case {'mean', 'bin'}
            if length(procpar.datadef(i).par) >= 1
                n = length(procpar.datadef(i).par{2});
            else
                n = nmeas;
            end
            scan.loops(loop).procfn(i+npc).dim = [n, procpar.n(2)*nrep];

        case 'hist'
            if length(procpar.datadef(i).par) >= 2
                n = length(procpar.datadef(i).par{2});
            else
                n = nmeas;
            end
            scan.loops(loop).procfn(i+npc).dim = [n, size(procpar.datadef(i).par{1}, 2)-1];
            
        case {'mom', 'bmom'}
            if isempty(procpar.datadef(i).par{1})
                scan.loops(loop).procfn(i+npc).dim = nmeas;
            else
                scan.loops(loop).procfn(i+npc).dim = size(procpar.datadef(i).par{1}, 2);
            end

        case {'corr', 'bcorr'}
            if isempty(procpar.datadef(i).par{1})
                n = nmeas;
            else
                n = size(procpar.datadef(i).par{1}, 2);
            end
            scan.loops(loop).procfn(i+npc).dim = [n, 2*procpar.datadef(i).par{2}+1];
            
    end
end