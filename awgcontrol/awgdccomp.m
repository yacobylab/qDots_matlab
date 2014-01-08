function awgdccomp(measp, workp, shiftscale, qpcfn, pgind, mplist, pplist, loop, inst)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


% wait command
% test with HW sync pulses

global awgdata;
global smdata;

if nargin < 5 || isempty(pgind)
    pgind = 1:length(awgdata.pulsegroups);
end

if nargin < 6 
    mplist = [];
end

if nargin < 7 
    pplist = [];
end

if nargin < 8
    loop = inf;
end

dacch = smchanlookup({'GateR', 'GateL', 'GateQPC'});
dacic = smchaninst(dacch);

rng = [smdata.inst(dacic(1, 1)).data.rng(ceil(dacic(1, 2)/2), :); ...
    smdata.inst(dacic(2, 1)).data.rng(ceil(dacic(2, 2)/2), :);
    smdata.inst(dacic(3, 1)).data.rng(ceil(dacic(3, 2)/2), :)];

rngfac = vertcat(smdata.channels(dacch).rangeramp);

if nargin < 4 
    qpcfn = [];
end

if isempty(qpcfn)
    dacic(3, :) = [];
end

if floor((dacic(1, 2)-1)/8) ~= floor((dacic(2, 2)-1)/8)
    fmt = sprintf('B%d;C%d;D%%d;', [floor((dacic(:, 2)-1)/8), floor(mod(dacic(:, 2)-1, 8)/2)]');
else
    fmt = sprintf('C%d;D%%d;', floor(mod(dacic(:, 2)-1, 8)/2)');
end

if nargin < 9
    inst = smdata.inst(dacic(1, 1)).data.inst;
elseif ischar(inst)
    inst = fopen(inst, 'w');
end

fprintcheck(inst, '{');

[awgdata.pulsegroups.script] = deal(0);

labcount = 1;
for j = pgind;
    pg = awgdata.pulsegroups(j);
    % could also load from file.
    
    if any(pg.nrep == Inf) || pg.npulse == 0;
        continue;
    end
    
    npulse = min(loop, pg.npulse);
    pulsedata = awgdata.pulsedata(pg.pulses(1:npulse));
    measpls = false(1, npulse);
    for i = 1:npulse
        measpls(i) = (pulsedata(i).pulsetab(1, end) > 2 || any(pg.pulses(i)==mplist))...
            && ~any(pg.pulses(i)==pplist);
    end
    
    %dt = zeros(1, pg.npulse);

    shift = awgdata.shift(:, pg.pulses(1:npulse))*shiftscale;
    val = zeros(size(shift)+[~isempty(qpcfn), 0]);

    val(1:2, :) = shift + repmat(workp', 1, npulse);
    val(1:2, measpls) = repmat(measp', 1, sum(measpls));% + shift(i, measpls(1));

    for i = 1:size(val, 2)*~isempty(qpcfn)
        val(3, i) = qpcfn(val(1:2, i)', smdata.chanvals);
    end
    
    for i = 1:size(val, 1)
        %shift(:, measpls) = shift(:, measpls(1));  
        val(i, :) = round((max(min(val(i, :), rngfac(i, 2)), rngfac(i, 1)) * rngfac(i, 4) - rng(i, 1))...
            ./diff(rng(i, :)) * 65535);
    end
    val = max(min(val, 65535), 0);

    
    fprintcheck(inst, '*%d:X%d;', labcount, labcount+2816);
    fprintcheck(inst, 'B%d;', floor((dacic(1, 2)-1)/8)); % redundant if different boards

    % need something here to wait for trigger if not provided by first pulse
    %fprintcheck(inst, '*%d:X%d;\n', labcount+1); % wait trig1 falling
    fprintcheck(inst, '*%d:', labcount+1);
    % this waits for a low on trig2
    
    awgdata.pulsegroups(j).script = labcount;
    labcount = labcount+2;
    
    for i = 1:npulse
        %hwsync = i < pg.npulse && size(pulsedata(i+1).marktab, 2) >= 1 && pulsedata(i+1).marktab(1, 1) == 0 ...
        %    && pulsedata(i+1).marktab(3, 1) ==  pulsedata(i+1).pulsetab(1, end);

        if i==npulse
            nextp = [i 1];
        else
            nextp = i+[0 1];
        end
        
        mark = [0 0];
        for k = 1:2
            if isempty(pulsedata(nextp(k)).marktab)
                mark(k) = 0;
            elseif pulsedata(nextp(k)).marktab(3, 1) == pulsedata(nextp(k)).pulsetab(1, end)
                mark(k) = 1;
            else
                mark(k) = -1;
            end
        end

        if any(mark == -1) || mark(1) == mark(2)
            fprintcheck(inst, '$%d;', round(pulsedata(i).pulsetab(1, end) * pg.nrep(min(i, end)) ...
                * 1e6 * pulsedata(i).tbase/pulsedata(i).clk));
            mark = [0, 0];
        end
        %if i == 1 || max(abs(shift(:, i) - shift(:, i-1))) > 1e-5
        % may be useful to improve timing or save memory
        fprintcheck(inst, fmt, val(:, i));
        %else accumulate time
        %end

        fprintcheck(inst, '*%d:X%d;', labcount, labcount+256*(3+mark*[8; 7]));
        labcount = labcount+1;
    end

    %fprintcheck(inst, 'X%d;', j*~isempty(strfind(pg.ctrl, 'loop')));
    fprintcheck(inst, 'X%d;', awgdata.pulsegroups(j).script+1);
end

fprintcheck(inst, '}');
if isnumeric(inst) && inst > 2
    fclose(inst);
end

if labcount > 256
    error('Out of DAC labels. Script invalid.');
end

function fprintcheck(inst, fmt, varargin)

str = sprintf(fmt, varargin{:});
fprintf(inst, '%s\n', str);
if ~isnumeric(inst) && ~strcmp(fscanf(inst), [str, 10]);
    error('Error writing script. Try again, flush buffer, look for a bug.');
end
