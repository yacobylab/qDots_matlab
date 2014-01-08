function awgmakegrp(pulsegroups, name, chans)
% awgmakegrp(pulsegroups, name)
% pulsegroups is a struct array with fields
% pulses and nrep, ctrl (optional);
% ctrl: loop (only for first loop)
% jump; array with from/to in row 1,2
% name (optional)
% chan ; number like like 34  or 21 to permute channels
% pulses and chan may be an array, in which case pulses in each row are
%   seperately permuted and summed; ie,
%   pulses = [ 7 8 ; 11 12 ]; chan = [ 12 ; 34 ] will apply pulse 7 the 8
%   to channel 1,2 while simultaneously running 11 then 12 on 34. Pulse
%   times must match up.
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

pls = [];
files={};
dir=[awgdata.datadir name];
if(fileattrib(dir) ~= 1)  
  mkdir(dir);
else
  if isfield(awgdata, 'grpdir')
    [oldgrp, oldchans]=load([awgdata.grpdir, name], 'pulsegroups', 'chans');
  else
    [oldgrp, oldchans]=save([awgdata.datadir, name], 'pulsegroups', 'chans');
  end
  if(~all(oldchans == chans))
      oldgrp=[];
      oldchans=[];
  end
end
trigname=sprintf('%s/trig',name);
offname=sprintf('%s/off',name);
awgmakepulse('write regen',awgdata.trigp,[],trigname);
awgmakepulse('write regen',awgdata.offp,[],offname);

zpulses=[];

zpulses=makezeros(zpulses,name,max(awgdata.pulsedata(awgdata.trigp).pulsetab(1,:)));
zpulses=makezeros(zpulses,name,max(awgdata.pulsedata(awgdata.offp).pulsetab(1,:)));

for i = 1:length(pulsegroups)
    pulsegroups(i).npulse = size(pulsegroups(i).pulses,2);
    pulsegroups(i).nchan = size(pulsegroups(i).pulses,1);
    pulsegroups(i).seqind = size(pls, 1)+1;
    mkdir(sprintf('%s%s/g%d',awgdata.datadir,name,i));
    zch=0;
    for j=1:pulsegroups(i).npulse
        pulsegroups(i).filenames{j}=sprintf('%s/g%d/p%d_%d',name,i,i,j);
        
        % Make a hybrid pulse, including all channel swapping.
        clear hb;
        for l=1:pulsegroups(i).nchan
            hb(l).pulses=pulsegroups(i).pulses(l,j);
            if isfield(pulsegroups(i),'chan') && ~isempty(pulsegroups(i).chan)
                if length(pulsegroups(i).chan) >= l
                    hb(l).chan = unpack_channels(pulsegroups(i).chan(l));
                else                    
                    hb(l).chan = unpack_channels(pulsegroups(i).chan(end));
                end
            else
              hb(l).chan=chans;            
            end
        end
        if(isfield(pulsegroups(i),'matrix'))
            matrix=pulsegroups(i).matrix;
        else
            matrix=[];
        end
        if(isfield(pulsegroups(i),'offset'))
            offset=pulsegroups(i).offset;
        else
            offset=[];
        end
        zcht=awgmakehybrid('write nostore',hb,pulsegroups(i).filenames{j},offset,matrix);
        pulsegroups(i).zch(j)=zcht;
        zch=bitor(zch,zcht);
    end
    fprintf('Pulsegroup %d (%15s); zch=%d\n',i,pulsegroups(i).name,zch);
    % Are there any zero channels?
    if(zch ~= sum(10.^(0:awgdata.nchan-1))) 
       tp = pulsegroups(i).pulses(1);
       l = max(awgdata.pulsedata(tp).pulsetab(1,:));
       zpulses=makezeros(zpulses,name,l);
    end
    % We now permute the pulse files with awgmakehybrid, not the columns
    % in the sequence file
    
    if isfinite(pulsegroups(i).nrep) % all finite
        if length(pulsegroups(i).nrep) == 1 && pulsegroups(i).npulse > 1
            if isfield(pulsegroups, 'jump') && ~isempty(pulsegroups(i).jump)
                pulsegroups(i).nrep = repmat(pulsegroups(i).nrep, 1, pulsegroups(i).npulse);
            else
                pulsegroups(i).nrep = [repmat(pulsegroups(i).nrep, 1, pulsegroups(i).npulse-1), 0];
                %Stay at last pulse unless group has jumps
            end
        end
        pulsegroups(i).npulse = size(pulsegroups(i).pulses(1,:),2);
        pls = [pls; [[awgdata.trigp; pulsegroups(i).pulses(1,:)'; awgdata.offp], [1;  pulsegroups(i).nrep'; 0], zeros(pulsegroups(i).npulse+2, 2), repmat(-1,pulsegroups(i).npulse+2,1), [11;pulsegroups(i).zch';11]]];
        files=[files {trigname}  pulsegroups(i).filenames {offname}];
        
    elseif length(pulsegroups(i).nrep) == 1
        pls = [pls; [pulsegroups(i).pulses(1,:)', zeros(pulsegroups(i).npulse, 3), repmat(-1,pulsegroups(i).npulse,1), pulsegroups(i).zch']];
        files=[files  pulsegroups(i).filenames ];
    else
        pulsegroups(i).nrep(~isfinite(pulsegroups(i).nrep)) = 0;
        pls = [pls; [pulsegroups(i).pulses(1,:)', pulsegroups(i).nrep', zeros(pulsegroups(i).npulse, 2), repmat(-1,pulsegroups(i).npulse,1), pulsegroups(i).zch']];
        files=[files pulsegroups(i).filenames];
    end
end

if(size(pls,1) > 4000)
    fprintf('Warning: too many lines (%d).  Aborting\n',size(pls,1));
    return;
else
    fprintf('%d lines in group\n',size(pls,1));
end

awgmakeseq('', pls, [name, '.seq'], chans,files, [name '/zero_%08d_1.wfm']);

if isfield(awgdata, 'grpdir')
    save([awgdata.grpdir, name], 'pulsegroups', 'chans');
else
    save([awgdata.datadir, name], 'pulsegroups', 'chans');
end
return;

function zpulses=makezeros(zpulses,name,l)
   % Make a new zero pulse.
    if(~any(abs(zpulses - l) < 1e-9))
      clear pinf;
      pinf.pulsetab = [0 l;0 0];
      pinf.marktab = [0; 0];
      pinf.tbase = 1000;
      pinf.name = sprintf('zero_%08d', round(pinf.tbase * pinf.pulsetab(1, 2)));
      awgmakepulse('write nostore', pinf, [], [name '/' pinf.name ]);
      zpulses(end+1)=l;
   end
return

function chans2=unpack_channels(c)
if(c ~= 0) % Cheap hack to pack channel list into 1 column
    chans2 = [];
    while (c ~= 0)
        chans2 = [mod(c,10) chans2];
        c=floor(c/10);
    end
else
    chans2 = [];
end
return

