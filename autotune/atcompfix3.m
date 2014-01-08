function grad=atcompfix3
%function grad=atcompfix3
% Fix the compensation matrix for one side or the other of the dot...
% requires an up-to-date charge scan.
% procedure: assume STP and TL peaks are a linear function of RF gates.
% accumulate STP and TL measurements into a 2xn_gate matrix ST
% ST and TL measurements with all gates zero is given by ST0.
% accumulate gate values into an n_gate x n_gate matrix., G
% then if (ST-ST0) = A . G,
% A, the reponse matrix, = (ST-ST0). G^(-1)
% once we have that, we can back out a corrected gradient for the response
% to the RF gates on the other dot.

% these should ideally be picked to not excite switches or charge changes
% in the dots.  order is RF gate order.

global tunedata;

if any(awgcntrl('israw')) 
  awgcntrl('amp');
end

gates={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'}; 

l = pdload('left'); r= pdload('right');
sides={'right','left'};
%inner gates have crazy pull with coupler, so reduce the voltage change;
grad=zeros(4);
for s=1:length(sides)
    smset(gates,[0 0 0 0]);
    dvl = -l.sep45.val*1e-3;
    dvr = -r.sep.val*1e-3;
    atswap(sides{s});
    if strcmp(sides{s},'left') %for side we run scans on, use small change
        dvl=[-.7e-3 .7e-3];
        dvr(1) = dvr(1); % go easy because of coupler
        dvr(2) = dvr(2);
    else
        dvr=[-.8e-3 .6e-3];
        dvl(2)=dvl(2)/2; dvl(1)=dvl(1)/2; 
    end
    % l.sep r.sep has form [2 1 3 4]
    curr_vals = cell2mat(smget(gates)); %find the gates current values -- all should be near 0. 
    curr_valsmat = repmat(curr_vals,4,1);
    scanvals=curr_valsmat'+[dvl(2) 0 0 0 ; 0 dvl(1) 0 0; 0 0 dvr(1) 0; 0 0 0 dvr(2)]';
    % accumulate responses in matrix of 2 x ngates, w/ row one for stp
    % response, row two for tl response. 
    st0(1,1)=stpscan(gates,curr_vals,sides{s},1);    st0(2,1)=tlscan(gates,curr_vals,sides{s},1);
    if(any(isnan(st0)))
        return;
    end
    for j=1:size(scanvals,2)
        stloc(1,j) = stpscan(gates,scanvals(:,j),sides{s},1+j);
        if(any(isnan(stloc)))
            return;
        end
        stloc(2,j) = tlscan(gates,scanvals(:,j),sides{s},1+j);
        if(any(isnan(stloc)))
            return;
        end
        stloc(:,j)=stloc(:,j)-st0;  % Offset from zero scan.
    end
    a=stloc*inv(scanvals-curr_valsmat'); % A is the response matrix.
    al=a(1:2,1:2); 
    ar=a(1:2,3:4);
    if strcmp(sides{s},'left');
        grad(1:2,1:2)=[-1 0 ; 0 -1];
        grad(3,1:2)=-pinv(al)*ar*[1 ; 0];
        grad(4,1:2)=-pinv(al)*ar*[0 ; 1];
    else
        grad(3:4,3:4)=[-1 0 ; 0 -1];
        grad(1,3:4)=-pinv(ar)*al*[1 ; 0];
        grad(2,3:4)=-pinv(ar)*al*[0 ; 1];       
    end
end
smset(gates,curr_vals);
grad=grad'
compmatrix_new=-pinv(grad)


load compmatrix compmatrix
fprintf('change in compmatrix \n');
compmatrix_new-compmatrix

info=input('(yes/no:', 's');
doit = strcmp(info, 'yes');
if doit    
    fname = ['compmatrix_' datestr(now,29)]; % see the help for why 29 
    save(fname,'compmatrix');
    compmatrix=compmatrix_new; 
    save compmatrix compmatrix; 
end
end

function stploc = stpscan2(gates,scanvals,side,ind)

global tunedata;
    nloop=1000;
    nrep = 5;

scan = fConfSeq2(tunedata.stp.plsgrp,{'nloop',nloop,'nrep',nrep, ...
    'datachan',tunedata.chrg.scan.loops(2).getchan,'opts','ampok'});

for i=1:length(gates)
    scan.consts(end+1).setchan=gates{i};
    scan.consts(end).val=scanvals(i);
end
data = smrun(scan); %, file);
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
subplot(3,2,ind);
pf=polyfit(x,d,1);
nd=smooth(d-pf(1)*x - pf(2));
ign=5;
[mm,mi]=max(nd(ign:end-ign));
p=fitwrap('plfit plinit samefig',x,d,[pf(2) mm x(mi+ign) range(x)/8 pf(1)],@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x,[0 1 1 1 0]);
p=fitwrap('plfit plinit samefig',x,d,p,@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x);
title(sprintf('ST+: %g; width %g',p(3)*1000,p(4)*1000));
stploc=p(3)*1000;
fprintf('STP loc is %g\n',stploc);
end
function tlloc = tlscan2(gates,gvals,side,ind)
global tunedata;
%file = smnext(sprintf('%s_tl_%02i',side,ind));
    nloop=1000;
    nrep=5;

scan = fConfSeq2(tunedata.tl.plsgrp,{'nloop',nloop,'nrep',nrep, ...
    'datachan',tunedata.chrg.scan.loops(2).getchan,'opts','ampok'});
for i=1:length(gates)
    scan.consts(end+1).setchan=gates{i};
    scan.consts(end).val=gvals(i);
end
data = smrun(scan); %, file);
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
subplot(3,2,ind);
p=fitwrap('plinit plfit samefig',x,d,[min(d), range(d), xm, range(x)/8, range(d)/2, xm+.1], ...
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
