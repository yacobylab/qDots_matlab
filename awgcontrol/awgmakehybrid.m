function zchan = awgmakehybrid(cntrl, pulses,filename_in, offset, matrix)
% zchan = awgmakehybrid(cntrl, pulses,filename, offset, matrix)
% 
%  Generates files for hybridized pulses.  Pulses is a struct array with
%  members:
%    pulses - pulse number
%    chan - channel map (ie [3 4]).  empty for identity.
%  Pulses are summed; ie including a pulse twice will double it's amplitude
% cntrl is a string controlling the behavior.
% reutrns a number indicating which channels are zero, ie 1100 means c1, c2
%   are zeros (read right-to-left)
% plot:     plot pulse
% write:    write the pulse  (filename must be set.)
%    offset is a row vector offset for the pulse.  Def 0
%    matrix is a rotation in  channel space to be applied before offset.
%      Def I

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

zchan=0;
global awgdata;

if ~isfield(awgdata, 'pulsedata')
    load(awgdata.datafile);
else
    pulsedata = awgdata.pulsedata;
end

% First figure out how big our output is
for q=1:length(pulses)
    plsnum=pulses.pulses;
    chan=pulses.chan;
    pulseinf = pulsedata(plsnum);   
    pulsetab = pulseinf.pulsetab;
    ninchan = size(pulsetab, 1)-1;

    pulses(q).npoints = round(max(pulsetab(1, :)) * pulseinf.tbase);
    pulses(q).mchan = max([ninchan pulses.chan]);
end

% Find nochan (number of output channels)
if(exist('matrix','var') && ~isempty(matrix))
    nochan=size(matrix,1);  
else
    nochan = max([pulses.mchan]);
end
data = zeros(nochan,max([pulses.npoints])+1);
chans=1:nochan;
marker = zeros(nochan,size(data,2)-1, 'uint8'); 

for q = 1:length(pulses)
    plsnum = pulses(q).pulses;
    pulseinf = pulsedata(plsnum);   
    pulsetab = pulseinf.pulsetab;
    nchan = size(pulsetab, 1)-1;

    npoints = round(max(pulsetab(1, :)) * pulseinf.tbase);
    % use one less? Maybe just skip first?
    time = linspace(pulsetab(1, 1), pulsetab(1, end), npoints+1);

    if pulseinf.taurc ~= Inf
        error('Taurc is broken for hybrid pulses');
    end
            
    if isempty(pulses(q).chan)
        chan=1:nchan;
    else
        chan=pulses(q).chan;
    end
    for j = 1:nchan
        cmask = ones(1,npoints+1);
        for i = 2:size(pulsetab, 2)
            mask = time >= pulsetab(1, i-1)-1e-11 & time <= pulsetab(1, i)+1e-11;
            % added small shifts to mitigate rounding errors 08/04/09. Never seen to matter.
            % below makes writes the pulse into data using lines to connect the
            % corners defined in pulstab
            mask = mask & cmask;
            data(chan(j), mask) = data(chan(j),mask)+(-pulsetab(j+1, i-1) * (time(mask) - pulsetab(1, i)) ...
                + pulsetab(j+1, i) * (time(mask) - pulsetab(1, i-1)))./...
                (pulsetab(1, i) -  pulsetab(1, i-1));
            cmask = cmask & ~mask;
        end
    end
    
    % lets pulses be defined with functions (eg. sin, cos) instead of
    % just lines
    for i = 1:length(pulseinf.pulsefn)
        mask = time > pulseinf.pulsefn(i).t(1) & time <= pulseinf.pulsefn(i).t(2);
        for j = 1:nchan
            data(chan(j), mask) = data(chan(j),mask)+pulseinf.pulsefn(i).fn{j}(time(mask)-pulseinf.pulsefn(i).t(1), pulseinf.pulsefn(i).args{j, :}) - avg(j);
        end
    end
        
    
    % extend marktab to be right dimensions but leave pulseinf.marktab intact
    % in order to write the original version in the database.
    % Also handle channel permutation.
    tempmarktab = pulseinf.marktab;
    tempmarktab(end+1:2*nchan+1,:) = 0; 
    
    for i = 1:size(tempmarktab, 2)
        for j = 1:nchan
            for k = 1:2;
                mask = time(1:end-1) >= tempmarktab(1, i) &...
                    time(1:end-1) < tempmarktab(1, i) + tempmarktab(j*2+k-1, i);
                marker(chan(j), mask) = bitor(marker(chan(j), mask), k);
            end
        end
    end
end

if exist('offset','var') && ~isempty(offset) 
  data=data+repmat(offset(:),1,size(data,2));
end

data(:, end) = [];

if exist('matrix','var') && ~isempty(matrix)
  data=matrix*data;
end


vc=data*0;
if any(abs(data(:)+vc(:)) > 1)
  fprintf('WARNING! Pulse exceeds range.\n');
end
    
if regexp(cntrl, 'plot')    
    figure(30);
    %clf;
    subplot(221)
    plot(time(1:end-1), data, time(1:end-1), round((data + vc)*2^13)./2^13);
    
    if nchan > 1
        subplot(222)
        plot(data(1, :), data(2, :), data(1, :) + vc(1, :), data(2, :) + vc(2, :));
    end
    
    if  any(marker(:))
        subplot(223)
        plot(time(1:end-1), [bitand(marker, 1); bitand(marker, 2)./2]);
    end
    return;
end
data = data + vc;
chans=0;
if ~isempty(strfind(cntrl, 'write'))
    for j = 1:nochan
        buf = reshape([reshape(typecast(single(data(j, :)), 'uint8'), 4, npoints); marker(j, :)],...
            5*npoints, 1);
        if(~all(buf == 0))
            zchan = zchan + 10^(j-1);
            filename = sprintf('%s_%d.wfm',filename_in,j);
            if 1 || strfind(cntrl, 'file')
                out = fopen([awgdata.datadir filename], 'wb', 'ieee-le');
                fwrite(out, [sprintf('MAGIC 1000\r\n#9%09d', 5*npoints), buf', ...
                    sprintf('CLOCK %.8e \r\n', pulseinf.clk)]);
                fclose(out);
            else
                %fwrite(out, [sprintf('MMEM:DATA "/wfm/%s", #7%07d', filename, 5*npoints + 45),...
                fwrite(out, [sprintf('MMEM:DATA "%s", #7%07d', filename, 5*npoints + 45),...
                    sprintf('MAGIC 1000\r\n#7%07d', 5*npoints), buf', sprintf('CLOCK %.8e \r\n', pulseinf.clk)]);
            end
        else
            %pulseinf.name
            %keyboard
        end
    end
end
if isempty(strfind(cntrl, 'regen')) && isempty(strfind(cntrl, 'nostore')) 
    if(~isempty(pulsedata))
      pulsedata(plsnum) = orderfields(pulseinf, pulsedata);    
    else
      pulsedata=[pulseinf];
    end
    if ~isfield(awgdata, 'pulsedata')
        save(awgdata.datafile, 'pulsedata');
    else
        awgdata.pulsedata = pulsedata;
    end
end
