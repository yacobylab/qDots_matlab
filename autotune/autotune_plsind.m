function autotune(ctrl, varargin)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


% autotune(ctrl, varargin)
% ctrl: 
%   ana: Analyse only, do not take data. varagin{1} = ind = run indecx
%   new: start new tuning series, varagin{1} = directory name
%   set: set new gate values as given in tunedata.runs(ind).gates
%   shift: compute new scan center based on previous run
%   all
%   chrg
%   lead
%   zoom
%   line
%   load 
%   read
%   grad: update derivative estimates.
%   plscal: oprional 3rd and 4th args are save index and list of
%            gates to measure (1 and 2).

global tunedata;
global smdata;
global awgdata;

% retake flag?

doscan = isempty(strfind(ctrl, 'ana'));

if strfind(ctrl, 'new')   
    tunedata.dir = varargin{1};
    mkdir(tuneedata.dir);
    tunedata.runs = [];
    ind = 1;
elseif doscan && length(varargin) < 1 
    ind = length(tundeata.runs) + 1;
else
    ind = varargin{1};    
end

doall = ~isempty(strfind(ctrl, 'all'));

if doscan
    if strfind(ctrl, 'set')
        smset(tunedata.gatechan, tunedata.runs(ind).gates);
    else
        tunedata.runs(ind).gates = cell2mat(smget(tunedata.gatechan));
    end
end

if doall || ~isempty(strfind(ctrl, 'shift'))  && ind > 1
    tunedata.cntr = (tunedata.runs(ind).gates-tunedata.runs(ind-1).gates) * tunedata.runs(ind-1).grad(:, 1:2) ...
        + tunedata.runs(ind-1).chrg(3:4) + tunedata.chrg.offset(:, 1)';
end

if doall || ~isempty(strfind(ctrl, 'chrg'))
    file = sprintf('%s/sm_chrg_%03i', tunedata.dir, ind);
    if doscan
        smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]); 
        tunedata.chrg.scan = smscanpar(tunedata.chrg.scan, tunedata.cntr);
        data = smrun(tunedata.chrg.scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data');
    end
    [tunedata.runs(ind).chrg, tunedata.runs(ind).slp] ...
        = chargeana(squeeze(mean(data{1}, 1)), vertcat(tunedata.chrg.scan.loops.rng));
 
    while 1
        str = input('Accept? (y/n)', 's');
        switch str 
            case {'y', 'Y'}
                break
            case {'n' 'N'}
                return
        end
    end
end

qpcfn = tunedata.chrg.scan.loops(2).trafofn{2};

if doall || ~isempty(strfind(ctrl, 'lead'))        
    for i = 1:length(tunedata.lead.pos)
        file = sprintf('%s/sm_lead%i_%03i', tunedata.dir, i, ind);
        if doscan
            if tunedata.lead.pos(i) > 0
                v = [-tunedata.runs(ind).slp(2); -1];
            else
                v = [1; 1/tunedata.runs(ind).slp(4)];
            end
            
            tunedata.measp(1:2, :) = tunedata.runs(end).chrg(3:4)' * [1, 1] + ...
                tunedata.lead.pos(i) * [v, v] + [[0; 0] tunedata.lead.bg(:, i)];
            for j = 1:size(tunedata.measp, 2)
                tunedata.measp(3, j) = qpcfn(tunedata.measp(1:2, j), smdata.chanvals);
            end
            tunedata.lead.scan.loops(2).trafofn = {@attrafo1, @attrafo2, @attrafo3};
            tunedata.lead.scan.loops(2).prefn(2).args = {awgseqind(-tunedata.lead.plsgrp(i))};

            data = smrun(tunedata.lead.scan, file);
            %plot(tunedata.measp(1, 1), tunedata.measp(2, 1), 'rx'); pause
            if any(isnan(data{1}(:))); return; end
        else
            load(file, 'data');
        end

        data = diff(squeeze(mean(data{1})));
        samprate = tunedata.lead.scan.consts(1).val;
        %nsamp = tunedata.lead.period * samprate *1e-6;
        x = (0:length(data)-1)./samprate * 1e6;

        tunedata.runs(ind).lead(i, :) = fitwrap('plinit plfit pause', x, ...tunedata.lead.period/nsamp * (0:nsamp-1), ...
            data, [mean(data), range(data)*sign(data(round(end/4))-data(round(3*end/4))), .3, .3, .1],  @leadfn);

    end
end

if doscan && (doall || ~isempty(strfind(ctrl, 'zoom')))     
    smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]);
    scan = smscanpar(tunedata.chrg.scan, tunedata.runs(end).chrg(3:4) + tunedata.chrg.offset,...
        tunedata.chrg.zoomrng, tunedata.chrg.zoomres);
    scan.loops(3:end) = [];
    %smset('PulseLine', awgseqind(tunedata.chrg.zoompls));
    scan.consts(2).val = awgseqind(tunedata.chrg.zoompls);
    data = smrun(scan, sprintf('%s/sm_zoom_%03i', tunedata.dir, ind));
    if any(isnan(data{1}(:))); return; end
end
  

if doscan && ( doall || ~isempty(strfind(ctrl, 'load')) || ~isempty(strfind(ctrl, 'read'))...
        || ~isempty(strfind(ctrl, 'line')))            
    tunedata.measp(1:2, :) = tunedata.runs(ind).chrg(3:4)' * [1, 1] + tunedata.offset;
    for j = 1:size(tunedata.measp, 2)
        tunedata.measp(3, j) = qpcfn(tunedata.measp(1:2, j), smdata.chanvals);
    end
    
    if doall || ~isempty(strfind(ctrl, 'zoom'))     
        figure(1000);
        subplot(122);
        hold on;
        plot(tunedata.measp(1), tunedata.measp(2), 'k.', 'markersize', 20);
        while 1
            str = input('Accept? (y/n)', 's');
            switch str
                case {'y', 'Y'}
                    break
                case {'n' 'N'}
                    return
            end
        end
    end
end

if doall || ~isempty(strfind(ctrl, 'line'))      
    file = sprintf('%s/sm_line_%03i', tunedata.dir, ind);
    if doscan        
        smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]);
        scan = smscanpar(tunedata.line.scan, tunedata.runs(end).chrg(3:4) + tunedata.chrg.offset);
            %tunedata.chrg.linerng, tunedata.chrg.lineres);
        %scan.loops(3:end) = [];
        %smset('PulseLine', awgseqind(awgdata.offp));
        scan.loops(2).trafofn{2} = qpcfn;
        data = smrun(scan, file);
        if any(isnan(data{1}(:))); return; end

        %tunedata.line.scan.loops(1).trafofn{2} = qpcfn;
        %tunedata.line.scan.trafofn = {@atlinetraf};
        %data = smrun(tunedata.line.scan, sprintf('%s/sm_line_%03i', tunedata.dir, ind));
    else
        load(file, 'data', 'scan');        
    end
    data = mean(data{1});
    x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
    tunedata.runs(ind).line = fitwrap('plinit plfit', x, data, [mean(data), 0, mean(x), .001 .0002], ...
        @(p, x)p(1)+p(2)*(x-p(3))-p(4)*tanh((x-p(3))./p(5)));
end

if doall || ~isempty(strfind(ctrl, 'load'))      
    file = sprintf('%s/sm_load_%03i', tunedata.dir, ind);
    if doscan
        tunedata.load.scan.loops(2).trafofn = {@attrafo1, @attrafo2, @attrafo3};
        data = smrun(tunedata.load.scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data');
    end
    data = -diff(squeeze(mean(data{1})));
    x = awgdata.xval(tunedata.load.scan.data.pulsegroups.pulses)*1e-3;

    tunedata.runs(ind).load = fitwrap('plinit plfit', x, ...
        data, [min(data), range(data), .05],  @(p, x)p(1)+p(2)*exp(-x./p(3)));
end

if doall || ~isempty(strfind(ctrl, 'read'))        
    file = sprintf('%s/sm_read_%03i', tunedata.dir, ind);
    if doscan
        tunedata.read.scan.loops(2).trafofn = {@attrafo1, @attrafo2, @attrafo3};
        data = smrun(tunedata.read.scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data');
    end

    data = -squeeze(diff(mean(data{1}), [], 3));
    data = data(1:end-1, :) - repmat(data(end, :), size(data, 1)-1, 1);
    samprate = tunedata.read.scan.consts(1).val;
    x = (0:length(data)-1)./samprate * 1e6;

    mask = round(tunedata.read.blank*samprate*1e-6):size(data, 2);

    tunedata.runs(ind).read = fitwrap('plinit plfit', x(mask), ...
        data(:, mask), [max(data(1, mask)), 10],  @(p, x)p(1)*exp(-x./p(2)));

end


if ~isempty(strfind(ctrl, 'grad'))  && ind > 1
    
    if isempty(tunedata.runs(ind).chrg)
        tunedata.runs(ind).vals(1:2) = nan;
    else
        tunedata.runs(ind).vals(1:2) = tunedata.runs(ind).chrg(3:4);
    end
    
    if isempty(tunedata.runs(ind).lead)
        tunedata.runs(ind).vals(3:4) = nan;
    else
        tunedata.runs(ind).vals(3:4) = mean(tunedata.runs(ind).lead(:, 3:4), 2)';
    end

    if isempty(tunedata.runs(ind).read)
        tunedata.runs(ind).vals(5) = nan;
    else
        tunedata.runs(ind).vals(5) = tunedata.runs(ind).read(2);
    end

    if isempty(tunedata.runs(ind).load)
        tunedata.runs(ind).vals(6) = nan;
    else
        tunedata.runs(ind).vals(6) = tunedata.runs(ind).load(3);
    end
    
    if isempty(tunedata.runs(ind).line)
        tunedata.runs(ind).vals(7) = nan;
    else
        tunedata.runs(ind).vals(7) = tunedata.runs(ind).line(5);
    end
    
    %tunedata.runs(ind).vals = [tunedata.runs(ind).chrg(3:4), mean(tunedata.runs(ind).lead(:, 3:4), 2)', ...
    %    tunedata.runs(ind).read(2), tunedata.runs(ind).load(3), tunedata.runs(ind).line(5)];
    
    % compute basis with one vector along last change, others orthogonal to it.
    d = tunedata.runs(ind).gates-tunedata.runs(ind-1).gates;
    n = sqrt(d*d');
    if n < 1e-3
        tunedata.runs(ind).grad = tunedata.runs(ind-1).grad;
    else
        d = d./n;

        M = eye(4) - diag(d) * repmat(d, 4, 1);
        [m, mind] = max(abs(d));
        M(mind, :) = d;

        % transform gradient to this basis
        grad = M * tunedata.runs(ind-1).grad;

        % set derivative along d
        grad(mind, :) = (tunedata.runs(ind).vals-tunedata.runs(ind-1).vals)/n;
        % transform back to old basis
        grad = inv(M) * grad;
        mask = isfinite(grad(:, 1));
        tunedata.runs(ind).grad(mask, :) = grad(mask, :);
    end    
end

if ~isempty(strfind(ctrl, 'plscal'))     
    if length(varargin) >=2 
        ind2 = varargin{2};
    else
        ind2 = ind;
    end
    if length(varargin) >=3 
        pos = varargin{3};
    else
        pos = 1:length(tunedata.plscal.pos);
    end
    for i = pos
        file = sprintf('%s/sm_plscal%i_%03i', tunedata.dir, i, ind2);

        if doscan
            smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]);

            if tunedata.plscal.pos(i) > 0 
                scan = smscanpar(tunedata.plscal.scan(i), ...
                    tunedata.runs(end).chrg(3:4) + tunedata.plscal.pos(i) * [-tunedata.runs(ind).slp(2), -1]);
                scan.loops(3-i).trafofn{2} = qpcfn;
            else
                scan = smscanpar(tunedata.plscal.scan(i), ...
                    tunedata.runs(end).chrg([4 3]) + tunedata.plscal.pos(i) * [1/tunedata.runs(ind).slp(4), 1]);
                scan.loops(3-i).trafofn{2} = @qpcfn2;
            end

            scan.consts(2).val = awgseqind(tunedata.plscal.pls(i));

            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end

        else
            load(file, 'data', 'scan');
        end
        data = mean(data{1});
        x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
        tunedata.plscal.fit(i, :) = fitwrap('plinit plfit', x, data, ...
            [mean(data), 0, mean(x), .001 .001 .002 .0002], ...
            @(p, x)p(1)+p(2)*(x-p(3))+p(4)*tanh((x-p(3)-p(6))./p(7))+p(5)*tanh((x-p(3)+p(6))./p(7)));

    end
    
    
end
end

function v = qpcfn2(x, y) % need this construct to avod saving too much workspace
global tunedata;
v = tunedata.chrg.scan.loops(2).trafofn{2}([0 x], y);
end
% function xit = exitcheck(data)
% if iscell(data)
%     data = data{1};
% end
% xit = any(isnan(data(:)));

function v = attrafo1(x, y)
global tunedata; v = tunedata.measp(1, x(2));
end

function v = attrafo2(x, y)
global tunedata; v = tunedata.measp(2, x(2));
end

function v = attrafo3(x, y)
global tunedata; v = tunedata.measp(3, x(2));
end

function v = atlinetraf(x, y)
global tunedata; v = tunedata.measp(1:2, 1)' + x(1) * [1 -1];
end

function y = leadfn(beta, x)
%global tunedata; 
%x0 = tunedata.lead.period/2;
x0 = x(end/2+1);
x = mod(x-beta(5), 2*x0);

y = beta(1) + .5 * beta(2) * ((cosh(.5 * x0/beta(3)) - exp((.5*x0-x)./beta(3)))./sinh(.5*x0/beta(3)) .* (x < x0) ...
    - (cosh(.5 * x0/beta(4)) - exp((1.5*x0-x)./beta(4)))./sinh(.5*x0/beta(4)) .* (x >= x0));
end
