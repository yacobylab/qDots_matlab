function [grad basis2]=atxyfix3(dv,side)
% function grad=atxyfix4(dv,side)
% Fix the gradient for the dot.
% Procedure: Assume STP and TL peaks are a linear function of gates.
%First, takes STP/TL scans to find the gradient for each qubit 
%(how much XL/YL/XR/YR move the triple point). This entails taking scans w/ no voltage
%then w/ voltage applied to each PlsRamp for that qubit. 
%Since chrg scans done w/ PlsRamps, a change in PlsRamp equivalent to a
%chng in triple point 
%Then find how XL, YL, XR, YR move STP/TL on that side, and normalize to how
%PlsRamps do. Do on both sides to get gradient.  
%With gradient, form basis usual way -- first invert, then mult. by current
%basis to get dependence on A gates. 
%stpscan and tl scan take a vector of gates and the voltages to apply to
%them -- in this case, all PlsRamps. 
%Note: based on how good the bases are, it one may want to change the
%vsteps, below, which decide how large of voltage to apply. In general,
%can't go too large or funny things will happen. 

%Note: there are many functions to create gradients. For the x,y gates, use atxyfix2 (for chrg
%scans), atxyfix3 (for stp/tl). For compensation matrix, use atcompfix3 (stp only). For other
%gates, use atgradfix2 (chrg scans), atgradfix3 (stp/tl). 

global tunedata;

if any(awgcntrl('israw')) 
  awgcntrl('amp');
end
gtvals=smget(1:19);

if ~exist('side','var')
    sides={'left','right'}; 
    gates={'XL','YL','XR','YR'}; 
    sing=0; 
else 
    sides={side};
    gates={};
    gates{1}=['X' upper(side(1))];
    gates{2}=['Y' upper(side(1))]; 
    sing=1; 
end
    b=basislookup(gates); 

rfs={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'}; 
if dv>2e-3
    dv=.2e-3; %The voltage change to be made on the PlsRamps
    fprintf('You have input a really large voltage. Reducing... \n')
end

figure(10); clf; figure(11); clf; figure(12); clf; figure(13); clf; 
grad=zeros(4); 

for s=1:length(sides)
   atswap(sides{s}); 
   stp0=stpscan(rfs,[0 0 0 0],sides{s},0); %No voltage applied on PlsRamps
   if isnan(stp0); smset(1:19,cell2mat(gtvals)); return; end
   tl0=tlscan(rfs,[0 0 0 0],sides{s},0);
   if isnan(tl0); smset(1:19,cell2mat(gtvals)); return; end
   
   if strcmp(sides{s},'left')
       vstep=[1 1 8 8].*dv; %apply smaller changes for gates closer to qubit
       chans=[2 1]; %the channels for PlsRamps on dot 
   else
       vstep=[6 3 1 1].*dv;
       chans=[3 4];
   end
   if sing
       vstep=[1 1].*dv;
   end
   for i=1:2
       rfv=zeros(1,4); rfv(chans(i))=dv; %apply dv to PlsRamp_chans(i). 
       comp(1,i) = (stpscan(rfs,rfv,sides{s},i)-stp0)/dv; % comp has format [stpX stpY; tlX tlY]
       comp(2,i) = (tlscan(rfs,rfv,sides{s},i)-tl0)/dv;
   end
   for j=1:length(gates) 
       atchg(gates{j},vstep(j));
       stppt=stpscan(rfs,[0 0 0 0],sides{s},j+2);
       tlpt=tlscan(rfs,[0 0 0 0],sides{s},j+2);
       atchg(gates{j},-vstep(j))
       dvec=[stppt-stp0 ; tlpt-tl0]; 
       grad(sort(chans),j)=-pinv(comp)*dvec./vstep(j);
       fprintf('A shift of %.1f mV on %s moves the STP and TL on %s by %3.1f uV and %3.1f uV, respectively \n',vstep(j)*1e3,gates{j},sides{s},dvec(1)*1e3,dvec(2));
   end
end
grad;
if sing
    grad=reshape(grad(find(grad)),2,2);
else
    chans=1:4;
end

basis=pinv(grad);
basis2=tunedata.basis(sort(chans),b)*basis;
printbasis(basis, basis2,gates,chans); 
info=input('(yes/no:', 's');
doit = strcmp(info, 'yes');
if doit
   tunedata.basis(sort(chans),b)=basis2;  
end

smset(1:19,cell2mat(gtvals))
close(10); close(11);  close(12); close(13); 
end


function printbasis(basis, basis2,oldbasis,chans)
    global smdata
    ngts=length(oldbasis); 
    newbasis={};
    for i=1:length(oldbasis)
        newbasis{i}=[oldbasis{i} 'n'];
    end
    basisprint=basis';
    basis2print=basis2';
    fprintf('       ')
    for i=1:ngts    
        fprintf(['%-4s',repmat('%-5s', 1, ngts+1), '\n'], '',oldbasis{i});
    end
    fprintf('\n');
    fprintf('-----------------------------------------------------\n')
    for i = 1:ngts
        fprintf(['%-9s:', repmat('%8.2g', 1, ngts), '\n'], newbasis{i}, basisprint(i,:));
    end      

    fprintf(['%-10s',repmat('%-8s', 1, ngts+1), '\n'], '',smdata.channels(chans).name);
    fprintf('\n');
    fprintf('-----------------------------------------------------\n')
    for i = 1:ngts
        fprintf(['%-9s:', repmat('%8.3g', 1, ngts), '\n'], newbasis{i}, basis2print(i,:));
    end
       
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