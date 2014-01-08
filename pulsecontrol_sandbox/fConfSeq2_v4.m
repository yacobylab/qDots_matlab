function scanSeq = fConfSeq2_v4(plsgrp, conf)
% function scanSeq = fConfSeq(plsgrp, conf)
% THIS VERSION ALLOWS GROUPS OF MULTI LENGTH PULSES
% Generate a scan using pulse groups plsgrp.
% conf can have any or none of:
%   nloop: default 200
%   oversamp: default 1
%   nrep: default 50.
%  datachan: default {'DAQ1','DAQ2'};
%   auxchan: default {'Time'};
%  fastmode: default 1
%   setmask:   default oversamp < 2
%  snglshot: default 2
%     0 - no snglshot  1 - smasnglshot   2 - built-in snglshot
% hwsampler: hardware sample rate, 40 Mhz default, nan to leave alone.
%      opts: random bool options.  Defaults to empty.
%            valid options are:
%              pol, fb, nodisp, nosave, ampok, nocheck,
%              raw, nosave
%            raw: save raw (non-averaged) data as well.
%            pol: add polairzed code
%             fb: add feedback control code
%         nodisp: disabel display
%          ampok: skip check that amplifier is off on awg
%        nocheck: skip sanity checks on groups
%         nosave: save only at end


if ~exist('conf','var')
    conf=struct();
end

if(iscell(conf))
    conf=struct(conf{:});
end
global awgdata;
global smdata;
global plsdata;
global fbdata;


conf.opts= def(conf,'opts','');  opts=conf.opts;
conf.nloop    = def(conf,'nloop',200); nloop = conf.nloop;
conf.oversamp = def(conf,'oversamp',1);  oversamp = conf.oversamp;
conf.nrep     = def(conf,'nrep',50);  nrep = conf.nrep;
conf.sampler  = def(conf,'sampler',nan);  sampler = conf.sampler;
conf.npulse   = def(conf,'npulse',nan);  npulse=conf.npulse;
conf.datachan = def(conf,'datachan',{'DAQ1','DAQ2'}); datachan = conf.datachan;
conf.auxchan  = def(conf,'auxchan',{'Time'}); auxchan = conf.auxchan;
conf.fastmode = def(conf,'fastmode',1);  fastmode=conf.fastmode;
conf.setmask  = def(conf,'setmask', oversamp < 2 ); setmask=conf.setmask;
conf.snglshot = def(conf,'snglshot',2); snglshot = conf.snglshot;
conf.hwsampler = def(conf,'hwsampler',40e6);  hwsampler = conf.hwsampler;
conf.extclk    = def(conf,'extclock',smdata.inst(14).data.extclk); extclk=conf.extclk;
if ischar(datachan)
    datachan = {datachan};
end
if ischar(auxchan)
    auxchan = {auxchan};
end

if isnan(oversamp) || oversamp==0
    oversamp = 1;
end

if feature('IsDebugMode')
   warning('debugger is on, fool!'); 
end

clear scanSeq;
scanSeq.data.conf=conf;  % Handy documentation.

plsgrp = awggrpind(plsgrp);

if plsgrp(1) > 0
    if isnan(npulse)
        npulse = awgdata(1).pulsegroups(plsgrp(1)).npulse(1);
    elseif npulse < 0
        npulse = awgdata(1).pulsegroups(plsgrp(1)).npulse(1) * abs(npulse);
    end

    if isempty(strfind(opts,'nocheck'))
      for ii = 1:length(plsgrp)
       if(npulse ~= awgdata(1).pulsegroups(plsgrp(ii)).npulse(1))
           warning('fConfSeq:PulseNum','Pulse number mismatch; length of group %s (%d) is not %d\n',...
               awgdata(1).pulsegroups(plsgrp(ii)).name,awgdata(1).pulsegroups(plsgrp(ii)).npulse(1),npulse);
       end
       if  plsinfo('stale',awgdata(1).pulsegroups(plsgrp(ii)).name)
           error('fConfSeq:StalePulse','Group needs uploading\n');
       end
%       plsmakegrp(plsgrp(ii),'check');  % Maybe we shouldn't do this -- it might be slow.
      end
    end
    
%     if isnan(sampler)
%         zl = plsinfo('zl', plsgrp(1));
%         sampler = awgdata(1).clk/(abs(zl(1, 1)) * max(1, awgdata(1).pulsegroups(plsgrp(1)).nrep(1)));       
%         
%         for ii=2:length(plsgrp)
%            zl = plsinfo('zl',plsgrp(ii));
%            sampler2 = awgdata(1).clk/(abs(zl(1, 1)) * max(1, awgdata(1).pulsegroups(plsgrp(ii)).nrep(1)));
%            if(sampler ~= sampler2)
%                warning('fConfSeq:PulseLen','Pulse length mismatch; sample rate for group %s is %g, not %g\n',...
%                    awgdata(1).pulsegroups(plsgrp(ii)).name,sampler2,sampler);
%            end
%         end
%     end    

    if isnan(sampler)
        zl = plsinfo('zl', plsgrp(1));
        % using max here is a hack to get the timing right for groups that
        % change pulse length. it means that we will acquire with the DAQ
        % longer than we want to . This should really get fixed to figure
        % it out properly
        sampler = awgdata(1).clk/(abs(max(abs(zl(:, 1)))) * max(1, awgdata(1).pulsegroups(plsgrp(1)).nrep(1)));       
        %sampler = awgdata(1).clk/(abs(mean(zl(:, 1))) * max(1, awgdata(1).pulsegroups(plsgrp(1)).nrep(1)));       
        
        for ii=2:length(plsgrp)
           zl = plsinfo('zl',plsgrp(ii));
           sampler2 = awgdata(1).clk/(abs(mean(zl(:, 1))) * max(1, awgdata(1).pulsegroups(plsgrp(ii)).nrep(1)));
           if(sampler ~= sampler2)
               warning('fConfSeq:PulseLen','Pulse length mismatch; sample rate for group %s is %g, not %g\n',...
                   awgdata(1).pulsegroups(plsgrp(ii)).name,sampler2,sampler);
           end
        end
    end    

    if (isfield(awgdata(1).pulsegroups(plsgrp(1)), 'jump') && ~isempty(awgdata(1).pulsegroups(plsgrp(1)).jump)) || ...
            awgdata(1).pulsegroups(plsgrp(1)).nrep(1) == 0 || ...  % single pulse repeated indefinitely -> same logic                   
            awgdata(1).pulsegroups(plsgrp(1)).nrep(1) == 1  % single pulse repeated indefinitely -> same logic                   
        npulse = npulse*nloop;
        fastmode = 1;
    else
        nloop = 0;
    end
    
    seqind = [awgdata(1).pulsegroups(plsgrp).seqind];
    scanSeq.data.pulsegroups = awgdata(1).pulsegroups(plsgrp);

else % single pulse, preceded by trigger.
    seqind = awgseqind(abs(plsgrp))-1; % subtract pulse to jump to preceding trigp.
    if isnan(npulse)
        npulse=1; 
    end  
    
    if isnan(sampler)
        pd = awgdata(1).pulsedata(abs(plsgrp(1)));
        sampler = pd.clk/(pd.pulsetab(1, end) * pd.tbase);
    end
    
    if nloop > 0 % treat same way as jump sequence
        npulse = npulse*nloop;
        fastmode = 1;
    end
    % could allow mix of groups and fixed pulses
end


scanSeq.loops(1).npoints = abs(npulse*oversamp);
scanSeq.loops(1).rng = [];

scanSeq.configfn.fn = @smabufconfig2;
%scanSeq.configfn.args = {'arm', 1:size(datachan)}; % used this line until 01/27/10. Doesn't make sense.
scanSeq.configfn.args = {'arm', 1}; % This configures only one channel. Seems to work.

scanSeq.loops(1).ramptime = 1/(oversamp * sampler);

scanSeq.loops(2).getchan = [datachan, auxchan];
scanSeq.loops(2).setchan = [];


scanSeq.loops(2).prefn(2).fn = '@(x, seqind, loop) smset(''PulseLine'', seqind(mod(floor(x(loop)-1), end)+1))';
scanSeq.loops(2).prefn(2).args = {seqind, 2};

if length(plsgrp) == 1 
    scanSeq.loops(2).npoints = nrep;
else
    scanSeq.loops(2).npoints = length(seqind);
    scanSeq.loops(3).npoints = nrep;
    scanSeq.loops(3).setchan = 'count';
end

if extclk
   if ~isnan(hwsampler)
     smset('RFfreq3',hwsampler);
     scanSeq.consts(1).setchan='RFfreq3';
     scanSeq.consts(1).val=hwsampler;
     scanSeq.consts(2).setchan='samprate';
     scanSeq.consts(2).val=hwsampler;
   end
   smdata.inst(14).data.extclk=1;
else
  if ~isnan(hwsampler)
     %smset('RFfreq3',10e6);
     smset('samprate',hwsampler);  % Yuck.  
     %scanSeq.consts(2).setchan='RFfreq3';
     %scanSeq.consts(2).val=10e6;
     scanSeq.consts(1).setchan='samprate';
     scanSeq.consts(1).val=hwsampler;   
  end
end
scanSeq = scanSeq.configfn.fn(scanSeq, scanSeq.configfn.args{:});

if scanSeq.loops(1).npoints ~= abs(npulse*oversamp)
    error('Invalid record length.');
    % checks for rounding error of machine (daq)
    % could be smarter here,  adapting numbers and/or cutting. Rate not
    % checked!
end

scanSeq.disp=[];
for l=1:length(datachan)    
  scanSeq.disp(end+1).loop = 2;
  scanSeq.disp(end).channel = l;
  scanSeq.disp(end).dim = 1;

  scanSeq.disp(end+1).loop = 2;
  scanSeq.disp(end).channel = l;
  scanSeq.disp(end).dim = 2;
end

if npulse*oversamp/nloop == 1
    scanSeq.disp(2) = []; % make disp work in degenerate case
end

switch fastmode
    case 1
        % remove inner loop for fast acq. 
        scanSeq.configfn.args{3} = [scanSeq.loops(1).npoints, 1/abs(scanSeq.loops(1).ramptime)];
        scanSeq.configfn.args{1} = 'fast arm';
        scanSeq.loops(1) = [];
        scanSeq.loops(1).prefn(2).args{2} = 1;
        scanSeq.loops(end).setchan = 'count';        
        if smdata.inst(sminstlookup('ATS660')).data.nrec == 0
            scanSeq.loops(1).waittime = -1;
            % we use async. reads now; no polling.  interrupt driven.            
        end        
        [scanSeq.disp.loop] = deal(1);
        
    case 2
        % shorten inner loop instead.
        scanSeq.configfn.args{3} = [scanSeq.loops(1).npoints, 1/abs(scanSeq.loops(1).ramptime)];
        scanSeq.configfn.args{1} = 'fast arm';
        scanSeq.configfn.args{4} = 2;
        scanSeq.loops(1).ramptime = scanSeq.loops(1).ramptime * scanSeq.loops(1).npoints;
        scanSeq.loops(1).npoints = 1;
        scanSeq.loops(end).setchan = 'count';
        scanSeq.saveloop = [length(scanSeq.loops) 5];
end
    

if nloop
    npulse = npulse/nloop; % restore previous value
    scanSeq.loops(1).procfn(1).fn(1).fn = @reshape;
    scanSeq.loops(1).procfn(1).fn.args = {[npulse * oversamp, nloop]};
    scanSeq.loops(1).procfn(1).fn(2).fn = @mean;
    scanSeq.loops(1).procfn(1).fn(2).args = {2};   
    scanSeq.loops(1).procfn(1).dim = npulse * oversamp;

    scanSeq.loops(1).procfn(2:length(datachan)) = scanSeq.loops(1).procfn(1);
    scanSeq.loops(1).procfn(length(datachan)+(1:length(auxchan))).fn = [];
end

if setmask && oversamp == 1 
    %pd = plsmakegrp(awgdata(1).pulsegroups(plsgrp(1)).name, 1);
    %&& isfield(pd.pulses(1).data, 'readout') && ~isempty(pd.pulses(1).data.readout)
    hwsr = cell2mat(smget('samprate'));
    rdw=[];
%     for i=1:length(plsgrp)
%         readout = plsinfo('ro', plsgrp(i));
%         for ll=1:2
%             l=find(readout(:,1) == ll);
%             if isempty(l)
%                 l = find(readout(:,1) == 0);
%             end
%             if isempty(l)
%                 l=1;
%             end
%             rdw(i,2*(ll-1)+(1:2))=readout(l(1),2:3);
%         end
%     end
rdout = plsinfo('ro',plsgrp(1));
for j = 1:length(plsgrp)
   rdout(j) = plsinfo('ro',plsgrp(j));
   masktmp = [];
   for k = 1:size(rdout(j).readout,3)
       if isfield(rdout(j),'plens')
        sdum = awgdata(1).clk/(abs(rdout(j).plens(k)) * max(1, awgdata(1).pulsegroups(plsgrp(1)).nrep(1)));
       else
          sdum = sampler;
       end
       dummy = false(2,hwsr/sdum);
       for ll = 1:2
         mind=round(hwsr * cumsum(rdout(j).readout(:,2:3,k))*plsdata.tbase/awgdata(1).clk);
         dummy(ll,mind(1):mind(2))=true;
       end
       masktmp = [masktmp, repmat(dummy,1,rdout(j).reps(k))];
   end
   masks{j} = (masktmp>0); %need to make it a logical
%    mind = [];
%    for k = 1:size(rdout(j).readout,3)
%        mind = [mind; round(hwsr * cumsum(rdout(j).readout(:,2:3,k))*plsdata.tbase/awgdata(1).clk)];
%    end
%    
 % maskind = round(hwsr * cumsum(readout(1:2))*plsdata.tbase/awgdata(1).clk);
end
%     [unique_masks, junk, masktab] =unique(rdw,'rows');
%     for i = 1:size(unique_masks,1)
%         masks{i} = false(2,hwsr/sampler);
%         for ll = 1:2
%           readout=unique_masks(i,2*(ll-1)+(1:2));
%           maskind = round(hwsr * cumsum(readout(1:2))*plsdata.tbase/awgdata(1).clk);
%           masks{i}(ll,maskind(1):maskind(2)) = 1;
%         end            
%     end
%     
    iic = smchaninst(datachan{1});
    scanSeq.loops(1).prefn(3:end+1)=scanSeq.loops(1).prefn(2:end);
    cntrlfn_str= func2str(smdata.inst(iic(1)).cntrlfn);
    scanSeq.loops(1).prefn(2).fn=['@(x,masks,loop) ',cntrlfn_str,'([smchaninst(''' datachan{1} '''),6],masks{mod(floor(x(loop)-1),end)+1})'];
    %scanSeq.loops(1).prefn(2).fn=['@(x,masks,loop) smcATS660v4_v2([smchaninst(''' datachan{1} '''),6],masks{mod(floor(x(loop)-1),end)+1})'];
    scanSeq.loops(1).prefn(2).args = {masks,1};
    if length(unique(rdout.plens))>1 % pulses changing lengths!
        fprintf('pulses changing lengths. this is a hack...\n');
        scanSeq.loops(1).prefn(3:end+1)=scanSeq.loops(1).prefn(2:end);
        % next line will set the num_pls_in_grp field in
        % smdata.inst(14).data so the daq driver knows what to do
        scanSeq.loops(1).prefn(2).fn = sprintf('@(x,t) %s([%d,%d,%d],t)',cntrlfn_str,iic(1),7,1); % 7 is this field 
        tmp = vertcat(awgdata.pulsegroups(plsgrp).npulse);
        scanSeq.loops(1).prefn(2).args = {tmp(:,1)};
    end
    %scanSeq.loops(1).prefn(2).fn=['@(x,masktab,masks,loop) smcATS660v4_v2([smchaninst(''' datachan{1} '''),6],masks{masktab(mod(floor(x(loop)-1),end)+1)})'];
    %scanSeq.loops(1).prefn(2).args={masktab, masks, 1};
    
    scanSeq.cleanupfn.fn = @smaconfigwrap;
    
    scanSeq.cleanupfn.args = {smdata.inst(iic(1)).cntrlfn, [iic, 6], []};    
    %scanSeq.cleanupfn.args = {@smcATS660v3, [smchaninst(datachan{1}), 6], []};    
end

%now set the n_pls_in_grp flag to empty
scanSeq.cleanupfn(end+1).fn = @smaconfigwrap;
scanSeq.cleanupfn(end).args = {smdata.inst(iic(1)).cntrlfn, [iic(1), 7, 1], []};   

if ~isempty(strfind(opts,'raw'))
    cpfn.fn = []; %copy function
    cpfn.args = {};
    ngc = length(scanSeq.loops(1).getchan);
    scanSeq.loops(1).procfn(1).fn = [cpfn, scanSeq.loops(1).procfn(1).fn];
    scanSeq.loops(1).procfn(1).fn(1).inchan = 1:length(datachan);
    scanSeq.loops(1).procfn(1).fn(1).outchan = length(scanSeq.loops(1).procfn) + (1:length(datachan));
    for k=scanSeq.loops(1).procfn(1).fn(1).outchan
        scanSeq.loops(1).procfn(k).fn=[];        
        scanSeq.loops(1).procfn(k).dim=nloop*npulse * oversamp;        
    end
 end

% time convention: time specifies end of sampling interval, i.e. time = (1:n)*dt,
% consistent with loadpulse.m. Want samples ending 

scanSeq.data.snglshot = snglshot;  % Later processing wants this.
switch snglshot 
    case 1 % use single shot framework to gather histogram and T1 info.
        clear pp;
        ndc = length(datachan);
        
        [pp.datadef(1:ndc).type] = deal('mom'); % std readout
        [pp.datadef(ndc + (1:ndc)).type] = deal('ghist'); %histograms

        for i = 1:ndc
            pp.datadef(i).par = {(1:npulse) + (i-1)*npulse};
            chan = smchaninst(datachan{i});
            pp.datadef(i+ndc).par = {linspace(-6e-3, 16e-3, 45) + fbdata.refval(chan(2)), (1:npulse) + (i-1)*npulse};
        end
    
   
        pp.datadef(2*ndc+1).type = 'ave'; % T1 measurement
        %pp.datadef(2*ndc+1).par = {'readout'};
        
        % this makes some of the above code obsolete, e.g. setting sample count in configfn and 
        % processing functions. Likely requires fastmode.
        scanSeq = smaSnglShot2(scanSeq, pp, 'raw dec', [], nloop);

        
    case 2 % Configure histogram "by hand"
        ngc = length(scanSeq.loops(1).getchan);
        scanSeq.loops(1).procfn(1).fn(2:end+1) = scanSeq.loops(1).procfn(1).fn(1:end);
        scanSeq.loops(1).procfn(1).fn(1).fn = [];
        scanSeq.loops(1).procfn(1).fn(1).args = {};        
        scanSeq.loops(1).procfn(1).fn(1).inchan = 1:length(datachan);
        scanSeq.loops(1).procfn(1).fn(1).outchan = length(scanSeq.loops(1).procfn) + (1:length(datachan));
        
        if(oversamp > 1)
            for i = 1:length(datachan)
                npts=4*45;
                scanSeq.loops(1).procfn(end+1).dim = [npts oversamp];
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @reshape;
                scanSeq.loops(1).procfn(end).fn(end).args = {[oversamp nloop]};
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @transpose;
                scanSeq.loops(1).procfn(end).fn(end).args = {};
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @histc;
                chan = smchaninst(datachan{i});
                % There is no compelling reason to limit the histogram to a
                % small region.  Disk space is cheap on this scale.
                %changing range around refval: 2013/02/26
                scanSeq.loops(1).procfn(end).fn(end).args = {linspace(-4*12e-3, 4*4e-3, npts) + fbdata.refval(chan(2))};
            end
        else
            for i = 1:length(datachan)
                npts=4*45;
                scanSeq.loops(1).procfn(end+1).dim = npts;
                scanSeq.loops(1).procfn(end).fn(1).fn = @histc;
                chan = smchaninst(datachan{i});
                % There is no compelling reason to limit the histogram to a
                % small region.  Disk space is cheap on this scale.
                %changing range around refval: 2013/02/26
                scanSeq.loops(1).procfn(end).fn(1).args = {linspace(-4*12e-3, 4*4e-3, npts) + fbdata.refval(chan(2))};
            end
        end
end 

if isempty(strfind(opts,'nosave'))
    scanSeq.saveloop = [length(scanSeq.loops), 5];
else
    scanSeq.saveloop = [-1, 1];
end
%check for channels in amp mode
if isempty(strfind(opts,'ampok'))
  israw = awgcntrl('israw');
  for i = 1:length(israw)
    if israw(i) ~=1
       fprintf('channel %d is in amp mode\n', i);
    end
  end
end

%check if AWGs on
if isempty(strfind(opts,'offok'))
  ison = awgcntrl('ison');
  for i = 1:length(ison)
    if ison(i) ~=1
       warning(sprintf('AWG %d is off\n', i));
       keyboard
    end
  end
end

if ~isempty(strfind(opts,'pol'))
   scanSeq = fScanPol(scanSeq);
elseif ~isempty(strfind(opts,'fb'))
   [t scanSeq] = fScanPol(scanSeq);
   fbdata.x=plsinfo('xval',plsgrp(1));
   fbdata.pulses=1:length(fbdata.x);
   %scanSeq.configfn(end+1).fn=@smaconfigwrap;
   %scanSeq.configfn(end).args{1}=@figure;
   %scanSeq.configfn(end).args{2}=1034;
   scanSeq.configfn(end+1).fn=@smaconfigwrap;
   scanSeq.configfn(end).args{1}=@feedback_control;   
   figure(1034);   
end

if ~isempty(strfind(opts,'predbz')) || ~isempty(strfind(opts,'logdbz'))
   a= struct('datachan',datachan);
   cc=struct('args',a,'name','pre_dbz','n_rep',1);
   if isfield(scanSeq,'configfn') && ~isempty(scanSeq.configfn) %want this to run before configuring daq card
       scanSeq.configfn(2:end+1) = scanSeq.configfn;
       scanSeq.configfn(1).fn = @sma_dBz_save;
       scanSeq.configfn(1).args = {cc};
   else %should really never happen
      scanSeq.configfn.fn = @sma_dBz_save;
      scanSeq.configfn.args = {cc};
   end
end

if ~isempty(strfind(opts,'postdbz')) || ~isempty(strfind(opts,'logdbz'))
   a= struct('datachan',datachan);
   cc=struct('args',a,'name','post_dbz','n_rep',1);
   scanSeq.cleanupfn(end+1).fn = @sma_dBz_save;
   scanSeq.cleanupfn(end).args = {cc};
end
% if isfield(conf,'fpga')
%     %perhaps silly and slow, but who cares. far easier than doing checking
%     conf.FPGA = conf.fpga;
% end
if ~isempty(strfind(opts,'FPGA')) || ~isempty(strfind(opts,'fpga'))
scanSeq.cleanupfn(end+1).fn = @get_FPGA_freqs;   
scanSeq.cleanupfn(end).args = {};
   if isempty(strfind(opts,'no_params'))
        scanSeq.configfn(end+1).fn = @store_FPGA_params;
        scanSeq.configfn(end).args = {};
   end
   scanSeq.loops(1).postfn(end+1).fn = @fpga_control_output;
   scanSeq.loops(1).postfn(end).args = {};
end

% if ~isempty(strfind(opts,'FPGA2'))
% %    scanSeq.cleanupfn(end+1).fn = @get_FPGA_data;
% %    scanSeq.cleanupfn(end).args = {};
%    if isempty(strfind(opts,'no_params'))
%        scanSeq.configfn(end+1).fn = @store_FPGA_params;
%        scanSeq.configfn(end).args = {};
%    end
% end


if ~isempty(strfind(opts,'nodisp'))
    scanSeq.disp=[];
end
return;

function v=def(s,f,v)
  if(isfield(s,f))
      v=getfield(s,f);
  end
return;

  