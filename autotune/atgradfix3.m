function atgradfix3(bases,dv,side)
%function atgradfix3(bases,dv,side)
%This is meant to improve a specific set of basis vectors, using the stp/tl format outlined in atcompfix3 and atxyfix3
%we normalize the changes using the changes of the X and Y gates 
% As in the the charge scan version, we then fix basis by adding some amount of
%X/Y bases to the gate basis that should cancel the movement of stp/tl
%scan. 
%Obviously, this relies on the X/Y bases being very good. 
%Starts by incrementing by small voltage, and as bases get better, doubles
%voltage change for each iteration. 
%bases can be 'all', 'leads', 'NTs', or any combination of gates, given in
%a script. 
%dv is the amoutn to change gates by. .1 to .4 e-3 is appropriate.
%if only one side is tuned up, can fix just one side. 

%Note: there are many functions to create gradients. For the x,y gates, use atxyfix2 (for chrg
%scans), atxyfix3 (for stp/tl). For compensation matrix, use atcompfix3 (stp only). For other
%gates, use atgradfix2 (chrg scans), atgradfix3 (stp/tl). 

global tunedata;
gtvals=smget(1:19);
[gates b]=getgates(bases); %puts selected bases into cell format. 
atswap('left'); atswap('right'); %to save tunedata
if ~exist('side','var')
    sides={'left','right'}; 
else 
    sides={side};
end
rfs={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'}; 
figure(10); clf; figure(11); clf; figure(12); clf; figure(13); clf; 

%For loop over the gates. For each gate, first run a charge scan to get original trip pts.
%Then increment by vstep, run scan for each side, and find how triple point
%moved. Update basis by adding in part of the X Y bases. 
for s=1:length(sides); 
    atswap(sides{s});
    stp0=stpscan(rfs,[0 0 0 0],sides{s},0); %No voltage applied on PlsRamps
    tl0=tlscan(rfs,[0 0 0 0],sides{s},0);
    xy{1}=['X' upper(sides{s}(1))]; 
    xy{2}=['Y' upper(sides{s}(1))]; 
    basxy=basislookup(xy); 
    for i=1:2
       atchg(xy{i},dv)
       comp(1,i) = (stpscan(rfs,[0 0 0 0],sides{s},i)-stp0)/dv; % comp has format [stpX stpY; tlX tlY]
       comp(2,i) = (tlscan(rfs,[0 0 0 0],sides{s},i)-tl0)/dv;
       atchg(xy{i},-dv)
    end
    for p=1:length(gates); 
        cont=1; vn=1; vstep=dv; 
        while cont
            atchg(gates{p},vstep);
            stppt=stpscan(rfs,[0 0 0 0],sides{s},p+2);
            if isnan(stppt); smset(1:19,cell2mat(gtvals)); return; end
            tlpt=tlscan(rfs,[0 0 0 0],sides{s},p+2);
            if isnan(tlpt); smset(1:19,cell2mat(gtvals)); return; end
            atchg(gates{p},-vstep)
            dvec=[stppt-stp0 ; tlpt-tl0]; 
            fprintf('A shift of %.1f mV on %s moves the STP and TL on %s by %3.1f uV and %3.1f uV, respectively \n',vstep*1e3,gates{p},sides{s},dvec(1)*1e3,dvec(2));
            grad=-pinv(comp)*dvec./vstep;
            dbasis=tunedata.basis(1:4,basxy)*grad;
            err(vn)=mean(grad);
            fprintf('Add %s %.2f, %s %.2f to %5s? \n',xy{1},grad(1),xy{2},grad(2),gates{p});   
            info=input('(yes/no:', 's');
            doit = strcmp(info, 'yes');
            if doit tunedata.basis(1:4,b(p))=tunedata.basis(1:4,b(p))+dbasis; end
            %the basis has format row 1:4 1a:4a, column 9:12 XL:YR
            vn=vn+1; 
            info=input('Try again? (yes/no): ','s');
            cont=strcmp(info,'yes'); 
            if cont 
                info=input('With a larger voltage? (yes/no): ','s');
                if strcmp(info,'yes')
                    vstep=2*vstep;           
                end
            end
        end
        fprintf('The error rate for %5s changed from %.2f to %.2f \n', gates{p},err(1),err(end)) 
    end
end

close(10); close(11); close(12); close(13); 
end

function [bases b]=getgates(opts)

bases=[];
if ~isempty(strfind(opts,'all'))
    bases={'Lead1', 'Lead2', 'Lead3', 'Lead4', 'T12', 'T34', 'N12', 'N34'};
end
if ~isempty(strfind(opts,'leads'))
    bases={'Lead1', 'Lead2', 'Lead3', 'Lead4'};
end
if ~isempty(strfind(opts,'NTs')) 
    bases= {'T12', 'T34', 'N12', 'N34'};
end

if ~isempty(strfind(opts,'Lead1')) bases{end+1}='Lead1'; end 
if ~isempty(strfind(opts,'Lead2')) bases{end+1}='Lead2';  end 
if ~isempty(strfind(opts,'Lead3')) bases{end+1}='Lead3';  end 
if ~isempty(strfind(opts,'Lead4')) bases{end+1}='Lead4';  end 
if ~isempty(strfind(opts,'T12')) bases{end+1}='T12';  end 
if ~isempty(strfind(opts,'T34')) bases{end+1}='T34';  end 
if ~isempty(strfind(opts,'N12')) bases{end+1}='N12';  end 
if ~isempty(strfind(opts,'N34')) bases{end+1}='N34';  end 

bases=unique(bases); 
b=basislookup(bases);
end        

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
subplot(4,6,ind+1);
pf=polyfit(x,d,1);
nd=smooth(d-pf(1)*x - pf(2));
ign=5;
[mm,mi]=max(nd(ign:end-ign));
p=fitwrap('plfit plinit samefig',x,d,[pf(2) mm x(mi+ign) range(x)/8 pf(1)],@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x,[0 1 1 1 0]);
p=fitwrap('plfit plinit samefig',x,d,p,@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x);
title(sprintf('ST+: %g; width %g',p(3),p(4)));
stploc=p(3);
fprintf('STP loc is %g\n',stploc*1000);
end

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
    subplot(4,6,ind+1);
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




function b=basislookup(basis)
%give this either a single gate in string form. 
global tunedata 
if ischar(basis)
    basis={basis}; 
end
if iscell(basis)
    for j = 1:length(basis); 
        b(j)=find(strncmp(basis{j},tunedata.basenames,length(basis{j}))); 
    end
else 
    fprintf('Please input a cell or a string')
    b=nan; 
end
end
