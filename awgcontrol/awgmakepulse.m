function plsnum = awgmakepulse(cntrl, pulseinf, plsnum,filename_in)
% plsnum = awgmakepulse(cntrl, pulseinf, plsnum,filename)
% 
% Adds awg pulses to pulsedata.   Pulses are not written to the instrument
% until they are added to a group.
% The parameters for a newly created pulse are stored in the struct
% array pulsedata in the file awgdata.datafile.
% 
% pulseinf is a struct parametrizing the pulse with the following fields:
% pulsetab: table specifying corner points of the pulse sequence.
%           The first row specifies the time in units of tbase clock cycles, the 
%           remainig ones the voltages in each channel. For points with the same 
%           time value, the voltage of the last point is used.
% tbase:    number of clock cycles per unit of time, default = 1000 (i.e. time in mus)
% clk:      clock rate (Hz), default 1e9.
% taurc:    RC time constant for bias filter compensation, default = Inf.
%           Inf suppresses average subtraction. Use -Inf to enable without correction.
% name:     user part of file name
% marktab:  matrix with marker information. First row contains switch on times, 
%           next two or 4 rows the pulse duration for CH1 and CH2 markers.
%           (this will need to change for N channel- might even need to
%           compensate for more than 2 markers per channel. also allow for
%           use of only some of the markers per channel)
%
% cntrl is a string controlling the behavior.
% regen     regenerate files for pulse number pulseinf from file
% awgdata.datafile
% nostore ; don't put in memory.
% plot:     plot pulse
% write:    write the pulse file.  (filename must be set.)
%
% The return value plsnum is the pulse index for the pulse.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if ~isfield(awgdata, 'pulsedata')
    load(awgdata.datafile);
else
    pulsedata = awgdata.pulsedata;
end

if regexp(cntrl, 'regen', 'once')
    plsnum = pulseinf;
    pulseinf = pulsedata(plsnum);   
    name = pulseinf.name;
else
    if ~isfield(pulseinf, 'name')
        pulseinf.name = '';
    end
        
    if ~isfield(pulseinf, 'tbase') || isempty(pulseinf.tbase)
        pulseinf.tbase = 1000;
    end

    if ~isfield(pulseinf, 'clk') || isempty(pulseinf.clk)
        pulseinf.clk = 1e9;
    end
    
    if ~isfield(pulseinf, 'taurc') || isempty(pulseinf.taurc)
        pulseinf.taurc = Inf;
    end

    if ~isfield(pulseinf, 'xval') || isempty(pulseinf.xval)
        pulseinf.xval = nan;
    end

    if ~isfield(pulseinf, 'marktab')
        pulseinf.marktab = [];
    end

    if ~isfield(pulseinf, 'pulsefn')
        pulseinf.pulsefn = [];        
    elseif ~isempty(pulseinf.pulsefn) && ~isfield(pulseinf.pulsefn, 'args')
        [pulseinf.pulsefn.args] = deal(cell(2, 0));
    end
   
     if ~isfield(pulseinf, 'readout')
        pulseinf.readout = [];
    end

    if nargin < 3
        plsnum = length(pulsedata) + 1;
    end
end

name = pulseinf.name;

if ~isempty(name)
    name = ['_', name];
end

if isfield(awgdata, 'plsgen');
    name = [awgdata.plsgen, name];
end

pulsetab = pulseinf.pulsetab;

nchan = size(pulsetab, 1)-1;

npoints = round(max(pulsetab(1, :)) * pulseinf.tbase);
% use one less? Maybe just skip first?

data = zeros(nchan, npoints+1);
time = linspace(pulsetab(1, 1), pulsetab(1, end), npoints+1);


if pulseinf.taurc == Inf
    avg = zeros(nchan, 1);
else
    avg = 0.5 * sum((pulsetab(2:end, 2:end) + pulsetab(2:end, 1:end-1)).*...
        repmat((pulsetab(1, 2:end) - pulsetab(1, 1:end-1)), nchan, 1), 2)./(pulsetab(1, end) - pulsetab(1, 1));
end

pulsetab(2:end, :) = pulsetab(2:end, :) - repmat(avg, 1, size(pulsetab, 2));
        
for j = 1:nchan
    for i = 2:size(pulsetab, 2)        
        
        mask = time >= pulsetab(1, i-1)-1e-11 & time <= pulsetab(1, i)+1e-11;
        % added small shifts to mitigate rounding errors 08/04/09. Never seen to matter.
        % below makes writes the pulse into data using lines to connect the
        % corners defined in pulstab
        data(j, mask) = (-pulsetab(j+1, i-1) * (time(mask) - pulsetab(1, i)) ...
            + pulsetab(j+1, i) * (time(mask) - pulsetab(1, i-1)))./...
            (pulsetab(1, i) -  pulsetab(1, i-1));
        
    end
end
        % lets pulses be defined with functions (eg. sin, cos) instead of
        % just lines
for i = 1:length(pulseinf.pulsefn)
    mask = time > pulseinf.pulsefn(i).t(1) & time <= pulseinf.pulsefn(i).t(2);
    for j = 1:nchan        
        data(j, mask) = pulseinf.pulsefn(i).fn{j}(time(mask)-pulseinf.pulsefn(i).t(1), pulseinf.pulsefn(i).args{j, :}) - avg(j);
    end
end
    
data(:, end) = [];
    %below calculates input voltage based on output voltage (different bc
    %of bias T.  see note(?)
if length(pulseinf.taurc) == 1
    vc = cumsum(data, 2) * (pulsetab(1, end) - pulsetab(1, 1))/(npoints * pulseinf.taurc);
else
    vc = cumsum(data, 2) * (pulsetab(1, end) - pulsetab(1, 1))./repmat(npoints * pulseinf.taurc', 1, npoints) ;
end

marker = zeros(nchan, npoints, 'uint8');

% extend marktab to be right dimensions but leave pulseinf.marktab intact
% in order to right the original version in the database
tempmarktab = pulseinf.marktab; 
tempmarktab(end+1:2*nchan+1,:) = 0; 

for i = 1:size(pulseinf.marktab, 2)        
    for j = 1:nchan
        for k = 1:2;
            mask = time(1:end-1) >= tempmarktab(1, i) &...
                time(1:end-1) < tempmarktab(1, i) + tempmarktab(j*2+k-1, i);
            marker(j, mask) = bitor(marker(j, mask), k);
        end
    end
end

if regexp(cntrl, 'plot')
    if any(abs(data(:)+vc(:)) > 1)
        fprintf('WARNING! Pulse exceeds range.\n');
    end
    
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

if ~isempty(strfind(cntrl, 'write'))
    for j = 1:nchan
        filename = sprintf('%s_%d.wfm',filename_in,j);
        buf = reshape([reshape(typecast(single(data(j, :)), 'uint8'), 4, npoints); marker(j, :)],...
            5*npoints, 1);

        if 1 || strfind(cntrl, 'file')
            out = fopen([awgdata.datadir filename], 'Wb', 'ieee-le');
            fwrite(out, [sprintf('MAGIC 1000\r\n#9%09d', 5*npoints), buf', ...
                sprintf('CLOCK %.8e \r\n', pulseinf.clk)]);
            fclose(out);
        else
            %fwrite(out, [sprintf('MMEM:DATA "/wfm/%s", #7%07d', filename, 5*npoints + 45),...
            fwrite(out, [sprintf('MMEM:DATA "%s", #7%07d', filename, 5*npoints + 45),...
                sprintf('MAGIC 1000\r\n#7%07d', 5*npoints), buf', sprintf('CLOCK %.8e \r\n', pulseinf.clk)]);
        end
    end
end

if isempty(strfind(cntrl, 'regen')) && isempty(strfind(cntrl,'nostore')) 
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
