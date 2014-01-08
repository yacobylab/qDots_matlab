function autotune(ctrl, varargin)

% autotune(ctrl, varargin)
% ctrl: 
%   ana: Analyse only, do not take data. varagin{1} = ind = run indecx
%   new: start new tuning series, varagin{1} = directory name
%   set: set new gate values as given in tunedata.runs(ind).gates
%   shift: compute new scan center based on previous run
%   all
%   chrg  ; x/y charge scan
%   ccopy ; copy previous charge scan.
%   lead  ; lead tunnel time by direct time domain charge sensing
%   zoom  ; locate measurement/load area
%   updateoff; change offset to match current measurement point
%   line  ; measure inter-dot coupling
%   load  ; measure load time by varying load pulse length
%   read  ; measure readout t1
%   stp: measure stp transition
%   tl: locate top lead
%   tmp: tune measurement point.
%   grad: update derivative estimates.
%   resp: Measure linear response of x,y transitions w.r.t. every gate.
%   copy grad: copy derivative estimates from prior run 
%   rgrad; override x,y portion of gradient matrix using results from resp.
%   plscal: oprional 3rd and 4th args are save index and list of
%            gates to measure (1 and 2).
%   
%   vargin is the run number.  If not specified, it defaults to last run
%     number, except for chrg and ccopy, which use the next run number.

global tunedata;
global smdata;
global awgdata;

% retake flag?

doscan = isempty(strfind(ctrl, 'ana'));

if strfind(ctrl, 'new')   
    tunedata.dir = varargin{1};
    mkdir(tunedata.dir);
    tunedata.runs = struct('gates', {}, 'chrg', {}, 'slp', {}, 'lead', {}, 'line', {}, ...
        'load', {}, 'read', {}, 'vals', {}, 'grad', {}, 'stp', {});
    ind = 1;
elseif length(varargin) < 1 || isempty(varargin{1}) 
    if (~isempty(strfind(ctrl,'chrg')) || ~isempty(strfind(ctrl,'ccopy')) || ~isempty(strfind(ctrl,'all'))) && isempty(strfind(ctrl,'ana'))
      ind = length(tunedata.runs) + 1;
      fprintf('Tune run %d %s starting\n',ind,tunedata.name);
    else
      ind = length(tunedata.runs);
      fprintf('Tune run %d %s continuing\n',ind,tunedata.name);
    end
else
    ind = varargin{1};    
end              

if ind > length(tunedata.runs)
    if(length(tunedata.runs) >= 1 && isempty(tunedata.runs(end).grad))
          q=input('Warning; gradient matrix not updated: Copy, Analyze, Ignore, Repeat ([c/A/i/r])? \n','s');
          if(isempty(q))
              q='a';
          end
          switch(q(1))
              case 'c'
                  autotune('grad copy ana',length(tunedata.runs));
              case 'r'
                  ind=length(tunedata.runs);
              case 'i'
              otherwise
                  autotune('grad ana',length(tunedata.runs));
%                  atprint('grad',length(tunedata.runs));
          end
    end
    tunedata.fine_ind=0;
end

doall = ~isempty(strfind(ctrl, 'all'));

if doscan
    if strfind(ctrl, 'set')
        smset(tunedata.gatechan, tunedata.runs(ind).gates);
    elseif ind > length(tunedata.runs)
        tunedata.runs(ind).gates = cell2mat(smget(tunedata.gatechan));
    elseif any(tunedata.runs(ind).gates ~= cell2mat(smget(tunedata.gatechan)));
        fprintf('WARNING: Gate values changed!\n');
    end
end


if (0 || ~isempty(strfind(ctrl, 'shift')))  && ind > 1
    if isempty(tunedata.chrg.offset)
        tunedata.cntr = (tunedata.runs(ind).gates-tunedata.runs(ind-1).gates) * tunedata.runs(ind-1).grad(:, 1:2) ...
            + mean(tunedata.runs(ind-1).chrg([1 2; 3 4])) ;
    else        
        tunedata.cntr = (tunedata.runs(ind).gates-tunedata.runs(ind-1).gates) * tunedata.runs(ind-1).grad(:, 1:2) ...
            + tunedata.runs(ind-1).chrg(3:4) + tunedata.chrg.offset(:, 1)';
    end
end

if ~isempty(strfind(ctrl,'ccopy'))
    tunedata.runs(ind).chrg = tunedata.runs(ind-1).chrg;
    tunedata.runs(ind).slp = tunedata.runs(ind-1).slp;    
end
%does charge scan and asks to click on the two triple pts
%chrgana2 analyzes the data to compute triple pts and slopes of transitions
runname = tunedata.name;
if ~isempty(runname)
    runname = ['_', runname];
end

if doall || ~isempty(strfind(ctrl, 'chrg'))
    file = sprintf('%s/sm_chrg%s_%03i', tunedata.dir, runname, ind);
    if doscan
        smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]); %un-obsolete
        scan = smscanpar(tunedata.chrg.scan, tunedata.cntr);
        data = smrun(scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data', 'scan');
    end
    if ndims(data{1}) == 3
        [tunedata.runs(ind).chrg, tunedata.runs(ind).slp] ...
            = chargean2(-squeeze(mean(data{tunedata.chrg.ch}, 1)), vertcat(scan.loops.rng));
    else
        [tunedata.runs(ind).chrg, tunedata.runs(ind).slp] ...
            = chargeana2(-data{tunedata.chrg.ch}, vertcat(scan.loops.rng));
    end
    if isempty(tunedata.runs(ind).chrg)
        return;
    end
    tunedata.measp(1:2, :) = tunedata.runs(ind).chrg(3:4)' * [1, 1] + tunedata.offset;
end

%copy the trafofn (qpcfn) into the correct loops of scan
if ~isempty(tunedata.chrg.scan.loops(1).trafofn)
    qpcfn = tunedata.chrg.scan.loops(1).trafofn{2};
elseif ~isempty(tunedata.chrg.scan.loops(2).trafofn)
    qpcfn = tunedata.chrg.scan.loops(2).trafofn{2};
else 
    qpcfn = @(x, y)0;
end

%takes a scan to compute tunneling times to leads
%pulses across a lead (gets measpt from location of TP & transition & slp)
% does same scan pulse away from lead to get bkgrnd measurement
if doall || ~isempty(strfind(ctrl, 'lead'))        
    for i = 1:length(tunedata.lead.pos)
        file = sprintf('%s/sm_lead%s_%i_%03i', tunedata.dir, runname, i, ind);
        if doscan
            mp = tunedata.measp;
            if tunedata.lead.pos(i) > 0
                v = [-tunedata.runs(ind).slp(2); -1];
            else
                v = [1; 1/tunedata.runs(ind).slp(4)];
            end
            
            tunedata.measp(1:2, :) = tunedata.runs(ind).chrg(3:4)' * [1, 1] + ...
                tunedata.lead.pos(i) * [v, v] + [[0; 0] tunedata.lead.bg(:, i)];
            for j = 1:size(tunedata.measp, 2)
                tunedata.measp(3, j) = qpcfn(tunedata.measp(1:2, j), smdata.chanvals);
            end
            tunedata.lead.scan.loops(2).trafofn = {@attrafo1, @attrafo2, @attrafo3};
            tunedata.lead.scan.loops(2).prefn(2).args = {awgseqind(-tunedata.lead.plsgrp(i))};
            tunedata.lead.scan.data.measp=tunedata.measp;
            data = smrun(tunedata.lead.scan, file);
            %plot(tunedata.measp(1, 1), tunedata.measp(2, 1), 'rx'); pause
            tunedata.measp = mp;
            if any(isnan(data{1}(:))); return; end
        else
            load(file, 'data');
        end

        data = diff(squeeze(mean(data{tunedata.lead.ch})));
        samprate = tunedata.lead.scan.consts(1).val;
        %nsamp = tunedata.lead.period * samprate *1e-6;
        x = (0:length(data)-1)./samprate * 1e6;
        try
            set(figure(501+i),'Name',sprintf('Lead %d',i));
            clf; hold on;
            tunedata.runs(ind).lead(i, :) = fitwrap('plinit plfit samefig', x, ...tunedata.lead.period/nsamp * (0:nsamp-1), ...
                data, [mean(data), range(data)*sign(data(round(end/4))-data(round(3*end/4))), .2, .2, .1],  @leadfn);
        catch
            fprintf('Fitting error\n');
        end
        switch(i)
            case 1
                l='x';
            case 2
                l='y';
            otherwise
                l=sprintf('%d',l);
        end
        fprintf('Lead %s: %g us, %g us\n',l,tunedata.runs(ind).lead(i,3),tunedata.runs(ind).lead(i,4));
        if strfind(ctrl, 'pause')
            pause;
        end
            
    end
end

% Crazy many-gated scan.  Generated by gradient_scans.m
if doall || ~isempty(strfind(ctrl, 'resp')) || ~isempty(strfind(ctrl,'cthulu'))
    filex = sprintf('%s/sm_resp%s_x_%03i',tunedata.dir, runname,ind);
    filey = sprintf('%s/sm_resp%s_y_%03i',tunedata.dir, runname,ind);
    s=0.9; % How far out to go in charge scan.
    if doscan
%        tunedata.resp.scanx.loops(1).rng=tunedata.chrg.scan.loops(1).rng;
        tunedata.resp.scanx.consts(2).val = tunedata.chrg.scan.loops(1).rng(1)*s;
        tunedata.resp.scany.consts(2).val = tunedata.chrg.scan.loops(2).rng(1)*s;
%        tunedata.resp.scany.loops(1).rng=tunedata.chrg.scan.loops(2).rng;
        
        datay=smrun(tunedata.resp.scany,filey); % run y-scan first. More likely to suck.
        datax=smrun(tunedata.resp.scanx,filex);
        if any(isnan(datax{1}(:))); return; end
        if any(isnan(datay{1}(:))); return; end
    else
        load(filex, 'data');
        datax=data;
        load(filey, 'data');        
        datay=data;
    end

    % Process and plot the data
    datax=squeeze(mean(datax{1},1));
    datay=squeeze(mean(datay{1},1));
    figure(77);
    clf;
    subplot(121);
    plot(datax');
    subplot(122);
    plot(datay');
    
    %Fit the data
    figure(78);
    clf;
    subplot(121);    
    tunedata.runs(ind).resp.gx=fit_resp(datax,tunedata.resp.scanx,ctrl);
    subplot(122);
    tunedata.runs(ind).resp.gy=fit_resp(datay,tunedata.resp.scany,ctrl);
    % Convert from our implicit dx, dy format to gradient matrix as expressed
    % by autotune.

    % dx_triple = -m2*dx/(m1-m2)-dy/(m1-m2)
    % dy_triple = m1*dx+dy

    %m2 = slope of nearly vert. line, m1=slope of nearly horz. line    
    xgate=find(strcmp(tunedata.resp.scanx.loops(2).setchan,tunedata.resp.scanx.loops(1).setchan));
    ygate=find(strcmp(tunedata.resp.scanx.loops(2).setchan,tunedata.resp.scany.loops(1).setchan));
    gradx=tunedata.runs(ind).resp.gx; % Gradients sweeping x
    grady=tunedata.runs(ind).resp.gy; % Gradients sweeping y
    gradx(xgate).val=-1;
    grady(ygate).val=-1;
    m1=grady(xgate).val
    m2=1./gradx(ygate).val
    tunedata.runs(ind).resp.tpgradx=-m2/(m1-m2)*[gradx.val]-[grady.val]/(m1-m2);
    tunedata.runs(ind).resp.tpgrady=m1*tunedata.runs(ind).resp.tpgradx+[grady.val];        
end

if doall || ~isempty(strfind(ctrl, 'zoom'))         
    for l=1:length(tunedata.chrg.zoompls)
      file=sprintf('%s/sm_zoom%s_%03i_%d', tunedata.dir, runname, ind,l);
      if(doscan)
        smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]);
        scan = smscanpar(tunedata.chrg.scan, tunedata.runs(ind).chrg(3:4) + tunedata.offset(:, 1)',...
          tunedata.chrg.zoomrng, tunedata.chrg.zoomres);
        % suspect was using a different offset before
        scan.loops(3:end) = [];
        scan.consts(2).val = awgseqind(tunedata.chrg.zoompls(l));      
        data = smrun(scan, file);
        if any(isnan(data{1}(:))); return; end        
      else
        load(file,'data','scan');
      end
      da{l}=data{1};
    end

    set(figure(2),'Name','Zoom Analysis');
    clf;
    sgn=-1;
    dn=da{1};
    for l=2:length(da)
     dn = dn+sgn*da{l};
     sgn=-sgn;
    end
    imagesc(scan.loops(1).rng,scan.loops(2).rng,dn);
    xlabel(scan.loops(1).setchan{1});
    ylabel(scan.loops(2).setchan{1});
    set(gca,'ydir','normal');
    axis image;
    hold on;

    h = plot(tunedata.measp(1, 1), tunedata.measp(2, 1), 'k.', 'markersize', 20);
    
    while 1;  
        str=input('[a]uto,[m]anual,accept(y/n)? ','s'); 
        switch str            
            case {'a', 'A'}
                delete(h);
                tunedata.measp(1:2, :) = tunedata.runs(ind).chrg(3:4)' * [1, 1] + tunedata.offset;
                h = plot(tunedata.measp(1, 1), tunedata.measp(2, 1), 'k.', 'markersize', 20);
            case {'y', 'Y'}
                break
            case {'m', 'M'}
                delete(h);
                tunedata.measp(1:2, :) = ginput(1)'* [1, 1] ...
                    + diff(tunedata.offset, [], 2) * [0, 1];
                h = plot(tunedata.measp(1, 1), tunedata.measp(2, 1), 'k.', 'markersize', 20);
            case {'n' 'N'}
                return  
        end
    end
    figure(1);
    hold on;
    plot(tunedata.measp(1,1),tunedata.measp(2,1),'k.','markersize',20);

    %for j = 1:size(tunedata.measp, 2)
    %    tunedata.measp(3, j) = qpcfn(tunedata.measp(1:2, j), smdata.chanvals);
    %end    
end

if doall || ~isempty(strfind(ctrl, 'line'))      
    file = sprintf('%s/sm_line%s_%03i', tunedata.dir, runname, ind);
    if doscan        
        smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]);
        cntr = mean(tunedata.runs(ind).chrg([1 2; 3 4]));
        % Permute this away from the center a bit to avoid the mousebite.
        cp = 0.33; % cp=0 means left triple, cp=1 means right.
        chrg = tunedata.runs(ind).chrg;
        cntr = [chrg(1)*cp+chrg(3)*(1-cp), chrg(2)*cp+chrg(4)*(1-cp)];
        %if isempty(tunedata.chrg.offset)
            %scan.loops(1).rng = scan.loops(1).rng + mean(tunedata.runs(ind).chrg([1 2; 3 4]));
        %    scan = smscanpar(tunedata.line.scan, cntr + [mean(tunedata.line.scan.loops(1).rng), 0]);
        %else
        %    scan = smscanpar(tunedata.line.scan, tunedata.runs(ind).chrg(3:4) + tunedata.chrg.offset);
        %end
        %scan.loops(2).trafofn{2} = qpcfn;
        
        %tunedata.chrg.linerng, tunedata.chrg.lineres);
        %scan.loops(3:end) = [];
        
        %for sweeping gate other than in charge scans, used for circular dot
        %scan = tunedata.line.scan;
        %vg = cell2mat(smget(scan.loops(1).setchan));
        %scan.loops(1).rng = scan.loops(1).rng + vg;
        %smset([tunedata.chrg.scan.loops(1).setchan(1:2), tunedata.chrg.scan.loops(2).setchan], ...
        %    [cntr(1), qpcfn(cntr), cntr(2)]);
        % assumes compensation in fast loop

        %data = smrun(scan, file);
        %smset(scan.loops(1).setchan, vg);
        %if any(isnan(data{1}(:))); return; end

        %tunedata.line.scan.loops(1).trafofn{2} = qpcfn;
        %tunedata.line.scan.trafofn = {@atlinetraf};
        %data = smrun(tunedata.line.scan, sprintf('%s/sm_line_%03i', tunedata.dir, ind));
        
        tunedata.line.scan.loops(1).trafofn(1).args = {cntr(1)};
        tunedata.line.scan.loops(1).trafofn(2).args = {cntr(2)};
        
        
        data = smrun(tunedata.line.scan, file); 
        if any(isnan(data{1}(:))); return; end        
        
    else
        load(file, 'data', 'scan');        
    end
    data = mean(data{tunedata.line.ch});
    x = linspace(tunedata.line.scan.loops(1).rng(1), tunedata.line.scan.loops(1).rng(2), tunedata.line.scan.loops(1).npoints);
    % do a linear fit for initial guess
    pts=1:(floor(length(x)/4));
    pf = polyfit(x(pts),data(pts),1);

    % Look at the skew of the derivative to guess a sign for the atan.
    dd = diff(data-pf(1)*x);
    skew = mean(dd)-median(dd);
    
    % Find the step
    [mm, mi] = max(abs(dd));
    cen=x(mi);
    sh=range(data-pf(1)*x)/2;
    %tunedata.runs(ind).line = fitwrap('plinit plfit woff', x, data, [pf(2)-pf(1)*(cen), pf(1), 0, .0002, .0002], ...
    %    @(p, x)p(1)+p(2)*(x-p(3))-p(4)*tanh((x-p(3))./p(5)));
    tunedata.runs(ind).line = fitwrap('plinit plfit woff fine', x, data, [pf(2)-pf(1)*(cen)+sh*sign(skew), +pf(1), cen, -sign(skew)*sh, range(x)/16.0], ...
        @(p, x)p(1)+p(2)*(x-p(3))-p(4)*tanh((x-p(3))./p(5)));
    fprintf('Line: %g meV\n',tunedata.runs(ind).line(5)*1e3);
end



if doall || ~isempty(strfind(ctrl, 'load'))      
    file = sprintf('%s/sm_load%s_%03i',tunedata.dir, runname, ind);
    if doscan
        %tunedata.load.scan.loops(2).trafofn = {@attrafo1, @attrafo2};
        tunedata.load.scan.consts(1).val = tunedata.measp(1,1);
        tunedata.load.scan.consts(2).val = tunedata.measp(2,1);
        tunedata.load.scan.data.pulsegroups = awgdata.pulsegroups(tunedata.load.plsgrp);
        tunedata.load.scan.loops(1).prefn(2).args{1} = awgseqind(-tunedata.load.plsgrp);

        data = smrun(tunedata.load.scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data');
    end
    if ndims(data{1}) == 3
        data = -diff(squeeze(mean(data{1})));
    else
        data = mean(data{1});
    end
    x = awgdata.xval(tunedata.load.scan.data.pulsegroups.pulses)*1e-3;

    tunedata.runs(ind).load = fitwrap('plinit plfit', x, ...
        data, [min(data), range(data), .05],  @(p, x)p(1)+p(2)*exp(-x./p(3)));
    fprintf('Load time: %g nsec\n',tunedata.runs(ind).load(3)*1e3);
end

if doall || ~isempty(strfind(ctrl, 'read'))        
    file = sprintf('%s/sm_read%s_%03i', tunedata.dir, runname, ind);
    if doscan
        tunedata.read.scan.consts(1).val = tunedata.measp(1,1);
        tunedata.read.scan.consts(2).val = tunedata.measp(2,1);
        tunedata.read.scan.data.pulsegroups = awgdata.pulsegroups(tunedata.read.plsgrp);
        tunedata.read.scan.loops(1).prefn(2).args{1} = awgseqind(-tunedata.read.plsgrp);
        data = smrun(tunedata.read.scan, file);
        if any(isnan(data{1}(:))); return; end
    else
        load(file, 'data');
    end

    set(figure(3),'Name','Readout T1');
    clf;
    subplot(121);
    plot(squeeze(mean(data{1}))');
        
    if ndims(data{1}) == 4
        data = -squeeze(diff(mean(data{1}), [], 3));
        data = data(1:end-1, :) - repmat(data(end, :), size(data, 1)-1, 1);    
    else
        set(figure(3),'Name','Readout T1');
        subplot(121);
        plot(squeeze(mean(data{1}))');
        data = -squeeze(diff(mean(data{1}), [], 2))';
    end
        
    samprate = tunedata.read.scan.consts(end).val;
    x = (0:length(data)-1)./samprate * 1e6;

    mask = round(tunedata.read.blank*samprate*1e-6):size(data, 2);

    dec=find(data(1,mask) < max(data(1,mask)*exp(-1)));
    if(~isempty(dec))
        lt=x(dec(1))-x(mask(1));
    else
        lt=10;
    end
    figure(3);
    subplot(122);    
    hold on;
    lt=3;
    tunedata.runs(ind).read = fitwrap('plinit plfit fine samefig', x(mask), ...
        data(:, mask), [max(data(1, mask))*1e3, lt],  @(p, x)1e-3*p(1)*exp(-x./p(2)));    
    fprintf('Readout T1: %g usec\n',tunedata.runs(ind).read(2));
end

if doall || ~isempty(strfind(ctrl,'stp'))
   file = sprintf('%s/sm_stp%s_%03i_%03i',tunedata.dir, runname,ind,tunedata.fine_ind);
   if(doscan)       
        tunedata.stp.scan.consts(1).val=tunedata.measp(1,1);
        tunedata.stp.scan.consts(2).val=tunedata.measp(2,1);
        tunedata.stp.scan.data.pulsegroups = awgdata.pulsegroups(tunedata.stp.plsgrp);
        tunedata.stp.scan.loops(1).prefn(2).args = {awgseqind(-tunedata.stp.plsgrp), 1};       
        
        % Hacking procfn for actual length.  Might be smarter to use
        % confSeq
        pgl=length(tunedata.stp.scan.data.pulsegroups.pulses);
        tunedata.stp.scan.loops(1).procfn(1).dim=pgl;
        tunedata.stp.scan.loops(1).procfn(1).fn(1).args{1}(1) = pgl;
        
        % Hack scan length and rate
        pd = awgdata.pulsedata(tunedata.stp.scan.data.pulsegroups.pulses(1));
        sampler = pd.clk/(pd.pulsetab(1, end) * pd.tbase * max(1, awgdata.pulsegroups(tunedata.stp.plsgrp).nrep(1)));       
        tunedata.stp.scan.configfn(1).args{3}(1:2)=[prod(tunedata.stp.scan.loops(1).procfn(1).fn(1).args{1}) sampler];
        
        data = smrun(tunedata.stp.scan, file);
        if any(isnan(data{1}(:))); return; end
   else
       load(file,'data');
   end    
   % Purely empirical fit form
   d=mean(data{1},1);
   x = awgdata.xval(tunedata.stp.scan.data.pulsegroups.pulses);
   set(figure(505),'Name','ST+ Analysis');
   clf;
   hold on;
   p=fitwrap('plfit plinit samefig',x,d,[median(d) range(d) mean(x) range(x)/2 ],@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)));
   tunedata.stp.eps=p(3);
   
   %tunedata.stp.stpdir is the direction and magnitude of step measured by xval
   tunedata.stp.stppnt=[tunedata.stp.scan.consts(1).val tunedata.stp.scan.consts(2).val] + p(3)*tunedata.stp.stpdir;
   tunedata.runs(ind).stppnt = tunedata.stp.stppnt;
   tunedata.runs(ind).stp_measp = tunedata.measp(1:2,1);
   tunedata.runs(ind).stp_eps = p(3);
   figure(2);
   hold on;
   plot(tunedata.stp.stppnt(1),tunedata.stp.stppnt(2),'kx','markersize',10);
   hold off;

   figure(1);
   hold on;
   plot(tunedata.stp.stppnt(1),tunedata.stp.stppnt(2),'kx','markersize',10);
   hold off;

end

if doall || ~isempty(strfind(ctrl,'tl'))
   file = sprintf('%s/sm_tl%s_%03i_%03i',tunedata.dir, runname, ind,tunedata.fine_ind);
   if(doscan)       
        tunedata.tl.scan.consts(1).val=tunedata.measp(1,1);
        tunedata.tl.scan.consts(2).val=tunedata.measp(2,1);
        tunedata.tl.scan.data.pulsegroups = awgdata.pulsegroups(tunedata.tl.plsgrp);
        tunedata.tl.scan.loops(1).prefn(2).args = {awgseqind(-tunedata.tl.plsgrp), 1};       
        
        % Hacking procfn for actual length.  Might be smarter to use
        % confSeq
        pgl=length(tunedata.tl.scan.data.pulsegroups.pulses);
        tunedata.tl.scan.loops(1).procfn(1).dim=pgl;
        tunedata.tl.scan.loops(1).procfn(1).fn(1).args{1}(1) = pgl;
        
        % Hack scan length and rate
        pd = awgdata.pulsedata(tunedata.tl.scan.data.pulsegroups.pulses(1));
        sampler = pd.clk/(pd.pulsetab(1, end) * pd.tbase * max(1, awgdata.pulsegroups(tunedata.tl.plsgrp).nrep(1)));       
        tunedata.tl.scan.configfn(1).args{3}(1:2)=[prod(tunedata.tl.scan.loops(1).procfn(1).fn(1).args{1}) sampler];
        tunedata.runs(ind).tl_measp = tunedata.measp(1:2,1);
        data = smrun(tunedata.tl.scan, file);
        if any(isnan(data{1}(:))); return; end
   else
       load(file,'data');
   end    
   % Purely empirical fit form
   x = awgdata.xval(tunedata.tl.scan.data.pulsegroups.pulses);
   d=mean(data{1},1);
   [cv ci]=max(d);
   xm=x(ci);
   set(figure(504),'Name','TR Analysis'); clf; hold on;
   p=fitwrap('plinit plfit samefig',x,d,[min(d), range(d), xm-range(x)/8, range(x)/8, range(d)/2, xm+range(x)/8], ...
       @(p,x) p(1)+p(2)*(tanh((x-p(3))/p(4))+1)/2 - p(5)*(tanh((x-p(6))/p(4))+1)/2);
   %tunedata.stp.stpdir is the direction and magnitude of step measured by
   %xval
   tl_de=(p(3)+p(6))/2;
   tunedata.runs(ind).tl_de=tl_de;
   fprintf('Relevant mismatch is %g\n',tl_de);  
end

% Fine-tune the measurement point using measured STP transititon and TR
% lead
if ~isempty(strfind(ctrl,'tmp'))
%     if(any(tunedata.measp(1:2,1) ~= tunedata.runs(ind).stp_measp) ||...
%             any(tunedata.measp(1:2,1) ~= tunedata.runs(ind).tl_measp))
%         error('Measurement point has changed. Will not work.\n');
%     end
   if any(tunedata.runs(ind).stp_measp ~= tunedata.runs(ind).tl_measp)
       error('Measurement different. Will not work.\n');
   end

  b1 = [-1 -tunedata.runs(ind).slp(tunedata.tmp.slp)]; b1 = b1 / norm(b1);
  b2 = tunedata.runs(ind).chrg(3:4)-tunedata.runs(ind).chrg(1:2); b2 = b2 / norm(b2);
  %bm = [ b1 ; b2 ]';
  %bmi = inv(bm);
  %bmi*b1';
  %bmi*b2';
  b1
  b2
  
  % Determine direction along b2 using lead data, along b1 using STP
  % stp_tgt is the target of the stp feedback
  c1 = (tunedata.runs(ind).stp_eps-tunedata.stp_tgt)*tunedata.stp.stpdir; % desired change
  c1 = c1 - dot(c1,b2)*b2;  % Distance perp to b2
  b1p = b1 - dot(b1,b2)*b2; % component of b1 perp to b2
  c1 = b1 * dot(c1,b1p) / (norm(b1p)^2);
  
  fprintf('Change for stp was %g,%g\n',c1(1),c1(2));
  
  c2 = (tunedata.runs(ind).tl_de-tunedata.tl_tgt)*tunedata.tl.tldir; % desired change
  c2 = c2 - dot(b1,c2)*b1;
  b2p = b2 - dot(b1,b2)*b1;
  c2 = b2 * dot(c2,b2p) / (norm(b2p)^2); 
  fprintf('Change for tl was %g,%g\n',c2(1),c2(2));

  big = norm(c1) + norm(c2) > 2e-3;
  while(big) % big change
     switch(input('Accept change (y/n)? ','s'))
         case 'y'
           big=0;
         case 'n'
           return;
     end             
  end
  tunedata.measp(1:2,1) = tunedata.runs(ind).tl_measp + c1' + c2';
  fprintf('New measurement point %f,%f\n',tunedata.measp(1,1),tunedata.measp(2,1));
  tunedata.fine_ind = tunedata.fine_ind + 1;
  figure(1);
  hold on;
  plot(tunedata.measp(1,1),tunedata.measp(2,1),'kv','markersize',10);
  
  figure(2);
  hold on;
  plot(tunedata.measp(1,1),tunedata.measp(2,1),'kv','markersize',10);

end

if ~isempty(strfind(ctrl,'updateoff'))
    tunedata.offset(1:2,1) = tunedata.measp(1:2,1)-tunedata.runs(ind).chrg(3:4)';
end

if isempty(tunedata.runs)
    return;
end

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
    
if ~isempty(strfind(ctrl, 'grad'))  && ind > 1
    %tunedata.runs(ind).vals = [tunedata.runs(ind).chrg(3:4), mean(tunedata.runs(ind).lead(:, 3:4), 2)', ...
    %    tunedata.runs(ind).read(2), tunedata.runs(ind).load(3), tunedata.runs(ind).line(5)];
    
    % Start by copying the old gradient.  That way, we'll keep "old" values
    % for any physical quant. not measured.
    tunedata.runs(ind).grad = tunedata.runs(ind-1).grad;
    % compute basis with one vector along last change, others orthogonal to it.
    d = tunedata.runs(ind).gates-tunedata.runs(ind-1).gates;
    n = sqrt(d*d');
    if n < .5e-3 || ~isempty(strfind(ctrl,'copy'))
      fprintf('Not recomputing gradient (%g < %g).\n',n,0.5e-3);
    else

        ng = length(tunedata.gradchan);

        if isfield(tunedata, 'basis')  % alternative basis given
            d2 = inv(tunedata.basis) * d';
            [m, mind] = max(abs(d2));
        end

        if isfield(tunedata, 'basis') && sum(abs(d2)) < m + 5e-4; 
            %change (approximately) along alternative basis vector.            
            M = tunedata.basis';
            n = d2(mind);
            fprintf('Using virtual gate basis.\n');
        elseif max(abs(d)) < .999 * n
            mind = [];
            fprintf('Copying gradient.\n');
            tunedata.runs(ind).grad = tunedata.runs(ind-1).grad;
        else
            d = d./n;
            M = eye(ng) - diag(d) * repmat(d, ng, 1);
            [m, mind] = max(abs(d));
            M(mind, :) = d;
            fprintf('Using physical gate basis.\n');
        end
        
        if ~isempty(mind)
            if ~isempty(tunedata.runs(ind-1).grad);
                % transform gradient to this basis
                grad(1:ng,:) = M(1:ng,1:ng) * tunedata.runs(ind-1).grad(1:ng,:);
            else
                grad = zeros(length(tunedata.gatechan), 7);
            end
            % set derivative along d
            grad(mind, :) = (tunedata.runs(ind).vals-tunedata.runs(ind-1).vals)/n;
            % transform back to old basis
            grad(1:ng,:) = inv(M(1:ng,1:ng)) * grad(1:ng,:);
            mask = isfinite(grad(1, :));
            tunedata.runs(ind).grad(1:ng, mask) = grad(1:ng, mask);
        end
    end  
    if(~isempty(strfind(ctrl,'rgrad')))
       for l=1:length(tunedata.gatechan)
           gn=smdata.channels(tunedata.gatechan(l)).name;
           gi=find(strcmp(tunedata.resp.scanx.loops(2).setchan,gn));
           if(isempty(gi))
               fprintf('Warning: gate %s not measured\n',gn);
           else
              tunedata.runs(ind).grad(l,1)=tunedata.runs(ind).resp.tpgradx(gi);
              tunedata.runs(ind).grad(l,2)=tunedata.runs(ind).resp.tpgrady(gi);       
           end
       end
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

            scan.consts(2).val = tunedata.plscal.pls(i);

            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end

        else
            load(file, 'data', 'scan');
        end
        data = mean(data{1});
        x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
        tunedata.plscal.fit(i, :) = fitwrap('plinit plfit', x, data, ...
            [mean(data), 0, mean(x), .001 .001 .0025 .0002], ...
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

% Fitting goodness
function [grad]=fit_resp(data,scan,ctrl)    
    fitfunc=@(p,x) p(1)+p(2)*x+p(3)*(1+tanh((x-p(4))/p(5)))/2;
    x=linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),scan.loops(1).npoints);    

    for i=1:size(data,1);
        ch=floor((i-1)/2)+1;
        fprintf('%d %s\n',i,scan.loops(2).setchan{ch});
        
        % offset slope step x0 width
        npt=15;
        dd=diff(data(i,:));       
        pf1=polyfit(x(1:npt),data(i,1:npt),1);
        pf2=polyfit(x(end-npt:end),data(i,end-npt:end),1);
        guess(2)=mean([pf1(1) pf2(1)]);
        guess(1)=mean(data(i,1:8)-guess(2)*x(1:8));
        ign=5;
        [m mi] = max(abs(smooth(dd(ign:end-ign))));
        mi=mi+ign;
        guess(3)=mean(data(i,end-8:end)-guess(2)*x(end-8:end) - guess(1));
        guess(4)=x(mi);
        guess(5)=range(x)/16;  
        guess=fitwrap('plinit plfit fine',x,data(i,:),guess, fitfunc);
        if(~isempty(strfind(ctrl,'pause')))
         pause;
        end
        figure(78);
        hold on;
        plot(x,ch+(data(i,:)-guess(1)-x*guess(2))/guess(3),'b+');
        plot(x,ch+(fitfunc(guess,x)-guess(1)-x*guess(2))/guess(3),'r-');
        if(mod(i,2) == 0)
            grad(ch).val=(guess(4)-ox)/scan.loops(2).trafofn(ch).args{2};
            grad(ch).name=scan.loops(2).setchan{ch};
        else
            ox=guess(4);
        end
    end    
end