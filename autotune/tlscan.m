function tlloc = tlscan(gates,gvals,side,ind)
%takes a tl scan w/ the offset given in scanvals, then fits to find tl
%pt.
    global tunedata;
    if isfield(tunedata.tl,'nloop')
        nloop = tunedata.tl.nloop;
    else
        nloop=200;
    end

    if isfield(tunedata.tl,'nrep')
        nrep = tunedata.tl.nrep;
    else
        nrep=5;
    end

    scan = fConfSeq2(tunedata.tl.plsgrp,{'nloop',nloop,'nrep',nrep, ...
        'datachan',tunedata.chrg.scan.loops(2).getchan,'opts','ampok'});
    for i=1:length(gates)
        scan.consts(end+1).setchan=gates{i};
        scan.consts(end).val=gvals(i);
    end
    data = smrun(scan);
    if any(isnan(data{1}(:))); tlloc=nan; return; end
    % Purely empirical fit form
    %x = awgdata.xval(tunedata.tl.scan.data.pulsegroups.pulses);
    x = plsinfo('xval', tunedata.tl.plsgrp);
    d=mean(data{1},1);
    [cv ci]=max(d);
    xm=x(ci);
    if strcmp(side,'left')
        figure(11); 
    else
        figure(13);
    end
    subplot(2,4,ind+1);
    p=fitwrap('plinit plfit samefig',x,d,[min(d), range(d), xm-range(x)/6, range(x)/8, range(d)/2, xm+range(x)/6], ...
        @(p,x) p(1)+p(2)*(tanh((x-p(3))/p(4))+1)/2 - p(5)*(tanh((x-p(6))/p(4))+1)/2);

    %tunedata.stp.stpdir is the direction and magnitude of step measured by

    %tl_de=(p(3)+p(6))/2;
    if p(6) > p(3)+ 400;
        tl_de = p(3)+200;
        fprintf('Broad peak. Using left edge + 200.\n')
    else
        tl_de = (p(3)+p(6))/2;
    end
    title(sprintf('TR %g',tl_de));
    fprintf('TL loc is %g\n',tl_de);
    tlloc=tl_de;
end