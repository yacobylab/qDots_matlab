function stploc = stpscan(gates,scanvals,side,ind)
%takes an stp scan w/ the offset given in scanvals, then fits to find stp
%pt. 
global tunedata;

if isfield(tunedata.stp,'nloop')
    nloop = tunedata.stp.nloop;
else
    nloop=1000;
end

if isfield(tunedata.stp,'nrep')
    nrep = tunedata.stp.nrep;
else
    nrep=5;
end

scan = fConfSeq2(tunedata.stp.plsgrp,{'nloop',nloop,'nrep',nrep, ...
    'datachan',tunedata.chrg.scan.loops(2).getchan,'opts','ampok'});

for i=1:length(gates)
    scan.consts(end+1).setchan=gates{i};
    scan.consts(end).val=scanvals(i);
end
data = smrun(scan);
if any(isnan(data{1}(:))); stploc=nan; return; end
% Purely empirical fit form
d = (mean(data{1},1));
%x = awgdata.xval(tunedata.stp.scan.data.pulsegroups.pulses);
x = plsinfo('xval', tunedata.stp.plsgrp);
if strcmp(side,'left')
    figure(10); 
else
     figure(12);
end
subplot(2,4,ind+1);
pf=polyfit(x,d,1);
nd=smooth(d-pf(1)*x - pf(2));
ign=5;
[mm,mi]=max(nd(ign:end-ign));
p=fitwrap('plfit plinit samefig',x,d,[pf(2) mm x(mi+ign) range(x)/8 pf(1)],@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x,[0 1 1 1 0]);
p=fitwrap('plfit plinit samefig',x,d,p,@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x);
stploc=p(3)*1e3;
stpw=p(4)*1e3;
title(sprintf('ST+: %g; width %g',stploc,stpw));

fprintf('STP loc is %g\n',stploc);
end