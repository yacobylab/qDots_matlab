function scan2 = qpccalib(scan, qrng, chans)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


% scan = qpccalib(scan, dind)
    
scan2 = scan;

% compensation loop is the one with trafofn's.
comploop = find(~cellfun(@isempty, {scan.loops.trafofn}), 1);
if isempty(comploop)
    error('No compensation gates found.')
end

if comploop == 2
    warning('comploop = 2');
end



gate1 = scan.loops(1).setchan;
if ~iscell(gate1)
    gate1 = {gate1};
end

gate2 = scan.loops(2).setchan;
if ~iscell(gate2)
    gate2 = {gate2};
end
ngates = [length(gate1), length(gate2)];
ngates(comploop) = sum(cellfun(@isempty, {scan.loops(comploop).trafofn.fn}));

gate1(ngates(1)+1:end) = [];
gate2(ngates(2)+1:end) = [];

gatecomp = scan.loops(comploop).setchan(ngates(comploop)+1:end);

nrep = 2;

if nargin < 2
    qrng = [.05, .001, 0];
end

if nargin < 3
    chans = 1:length(gatecomp);
end


scan.loops(2).setchan = {};
scan.loops(2).trafofn = {};
scan.loops(1).trafofn = {};
scan.loops(2).rng = [];
scan.loops(2).npoints = nrep;


plsscan = strcmp(gate2{1}(1:end-1), 'PlsDC');

if plsscan
    awgcntrl('exton', str2double(gate2{1}(end)));
    awgcntrl('extoff', 5-str2double(gate2{1}(end)));
    isrc = smquery('SR830', 'ISRC ?', '%s\n', '%d');
    exc = smget({'LockExc', 'LockSens'});

end


for dind = chans;

    vg = smget(gatecomp{dind});
    scan.loops(1).rng = vg{1} + qrng(1) * [-.5 .5] + qrng(3);
    scan.loops(1).setchan = gatecomp(dind);
    %scan.loops(1).ramptime = -.005;
    %scan.loops(1).npoints = 100;

    scan.cleanupfn.fn =  @smaconfigwrap;
    scan.cleanupfn.args = {@smset, gatecomp{dind}, [vg{:}]};

    vg = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
    %awgcntrl('exton', 2);
    scan.figure = 1001;
    d = smrun(scan);
    %awgcntrl('extoff', 2);
    scan.figure = 1002;
    if any(isnan(d{dind}(:)))
        return
    end
    
    
    d = mean(d{dind}, 1);
    %[m, mi] = max(mean(d{3}, 1)); % max slope - use if available
    %qpcc(2) = vg(mi);
   
    [m, mi] = max(d); % peak
    %mask = vg > vg(find(vg < vg(mi)-1e-3 & [1 diff(d) < 0], 1, 'last')) & vg < vg(mi);
    mask = find(vg < vg(mi)-5e-3 & [1 diff(d) < 0], 1, 'last'):mi;
    % positive slopes below peak
    %mask = vg > vg(mi)-.01 & vg < vg(mi);

    setp = mean([m, min(d(mask))]); % mean of min and max
    qpcc(1) = interp1(d(mask), vg(mask), setp);

    %dd = mean(d{dind});
    %mask = gradient(dd) < 0;
    %qpcc(1) = interp1(dd(mask), vg(mask), qpcdata.setp);
    if isnan(qpcc(1))
        return
    end

    mask = abs(vg-qpcc(1)) < 0.5 * qrng(2); 
    if plsscan && sum(mask) >= 5
        fqpc1 = d(mask)/[vg(mask); ones(1, sum(mask))];
        smset(gatecomp{dind}, qpcc(1));
    else
        % measure response to compensation gate over small range.
        scan.loops(1).rng = qpcc(1) + qrng(2) * [-.5 .5];
        scan.cleanupfn.args{3} = qpcc(1);
        d2 = smrun(scan);

        if any(isnan(d2{dind}(:)))
            return
        end

        fqpc1 = mean(d2{dind}/[linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);...
            ones(1, scan.loops(1).npoints)]);
    end
    
    
    % measure response to swept gate(s).
    if plsscan
        smprintf('SR830', 'ISRC 0');
        sl = cell2mat(smget('Lockin'));
        if dind == 2;
            smprintf('SR830', 'ISRC 1');
            sl = sl - cell2mat(smget('Lockin'));
        end
        smprintf('SR830', 'ISRC %d', isrc);
        awgcntrl('extoff',(str2double(gate2{1}(end))));

        qpcc(4) = -sl /(exc{1} * fqpc1(1)); % for reading directly

    else
        vg3 = smget(gate2);
        qpcc(5) = vg3{1};
        scan.loops(1).rng = vg3{1} + qrng(2) * [-1 1];
        scan.loops(1).setchan = gate2;
        scan.cleanupfn.args(2:3) = {gate2, [vg3{:}]};

        d3 = smrun(scan);

        if any(isnan(d3{dind}(:)))
            return
        end

        fqpc2 = mean(d3{dind}/[linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);...
            ones(1, scan.loops(1).npoints)]);
        qpcc(4) = -fqpc2(1)/fqpc1(1);
    end
    
    % fast gate
    vg4 = smget(gate1);
    qpcc(3) = vg4{1};
    scan.loops(1).rng = vg4{1} + qrng(2) * [-2 2];
    scan.loops(1).setchan = gate1;
    scan.cleanupfn.args(2:3) = {gate1, [vg4{:}]};

    d4 = smrun(scan);

    fqpc3 = mean(d4{dind}/[linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);...
        ones(1, scan.loops(1).npoints)]);

    %qqpc entries [y-offset, gate1 slope, gate1 x-offset, gate2slope, gate2 x-offset]
    % compensation fn = Yoffset +slope1*(x-Xoffset1)+slope2*(x-Xoffset2)
    qpcc(2) = -fqpc3(1)/fqpc1(1);
    scan2.loops(comploop).trafofn(ngates(comploop)+dind).fn = ...
        @(x,y,qpcc)qpcc(1)+qpcc(2)*(x(1)-qpcc(3)) + qpcc(4)*(x(2)-qpcc(5));
    scan2.loops(comploop).trafofn(ngates(comploop)+dind).args = {qpcc};
end
