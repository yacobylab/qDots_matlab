%% 
awgclk='AWGclock';
smset(awgclk,1.0e9);
smdata.inst(14).data.extclk = 1;

maxfreq=1e9; nfreq=200; drange=[1 3];
baseratio=10;
freqs=[maxfreq./linspace(drange(1),drange(2),nfreq)];


d=fConfSeq2([repmat(28,1,nfreq)],struct('nloop',50,'nrep',40,'opts','', 'datachan', 'DAQ2','hwsampler',freqs(1)/baseratio));
d.configfn(2)=d.configfn(1);
d.loops(2).setchan='count';
d.loops(2).npoints=50;
d.loops(1).setchan={'AWGclock','RFfreq3'};
d.loops(1).trafofn(1).fn=@(x,y) freqs(x(1));
d.loops(1).trafofn(1).args={};
d.loops(1).trafofn(2).fn=@(x,y) freqs(x(1))/baseratio;
d.loops(1).trafofn(2).args={};

if 0
  d.loops(1).postfn(2)=d.loops(1).postfn(1);
  d.loops(1).postfn(1).fn=@(x) smset('AWGclock',1e9);
  d.loops(1).postfn(1).args={};
end
d.saveloop=[2 1];

%Change the mask to be smart
mask=d.loops(1).prefn(2).args{2}{d.loops(1).prefn(2).args{1}(1)};
masks={mask};
if 0  % Only fix T1
    for j=2:length(freqs)
        masks{j}=mask*0;
    end
    for i=1:size(mask,1)
        first=find(mask(i,:),1,'first');
        len=find(mask(i,:),1,'last')-first;
        for j=2:length(freqs)
            masks{j}(i,first:(first+ceil(len*freqs(j)/1e9)))=1;
        end
    end
else %fix everythings
   ro=plsinfo('ro',d.data.pulsegroups(1).name);
   r=pdload('right');   
   for j=1:length(freqs)
       hwsr=freqs(j)/baseratio;
       masks{j}=false(size(mask));
       nro=ro;
       nro(2)=(nro(2)-r.meas.time(2))*freqs(1)/freqs(j); % Actual start of readout.
       nro(2)=nro(2)+r.meas.time(2);       
       maskind = round(hwsr*cumsum(nro(2:3))*plsdata.tbase/freqs(1));
       %fprintf('Frequency %g: Readout %g-%g (%d-%d)/%d\n',freqs(j),nro(2:3),maskind(1:2),size(masks{j},2));
       masks{j}(1:2,maskind(1):maskind(2))=1;
   end

end

%masktab
d.loops(1).prefn(2).args{1}=1:length(freqs);
%masks
d.loops(1).prefn(2).args{2}=masks;

d.data.frange = drange;
d.data.freqs = freqs;
d.cleanupfn(end+1).fn = @sleep;
d.cleanupfn(end).args = {};

'Scan ready'





%% A new variation; clock DAQ off Agilent generator. use FB
% Old file from Mikey do not run this part.
awgclk='AWGclock';
smset(awgclk,1.0e9);
alloff=awgseqind('all_off_LR');
dbzgrp = 28;
%exchgrp = 30; nfreq=800; drange=[1 3];
% time=[1 2 3 4 5 6 7 8 9 10]
%maxfreq=1e9; exchgrp = 24; nfreq=300; drange=[0.95/1 2.05/1];
%maxfreq=1e9; exchgrp = 25; nfreq=300; drange=[1.95/2 3.05/2];
%maxfreq=1e9; exchgrp = 26; nfreq=300; drange=[2.95/3 4.05/3];
maxfreq=1e9; exchgrp = 27; nfreq=300; drange=[3.95/4 3.05/4];


baseratio=10;
smdata.inst(14).data.extclk = 1;


freqs=[1e9 maxfreq./linspace(drange(1),drange(2),nfreq)];

d=fConfSeq2([dbzgrp repmat(exchgrp,1,nfreq)],struct('nloop',50,'nrep',50,'opts','pol', 'datachan', 'DAQ2','hwsampler',freqs(1)/baseratio));
%d=fConfSeq2([dbzgrp repmat(exchgrp,1,nfreq)],struct('nloop',100,'nrep',30,'opts','pol', 'datachan', 'DAQ2','hwsampler',100e6));
d.loops(1).setchan={'RFfreq3',awgclk};
d.loops(1).trafofn(2).fn=@(x,y) freqs(x(1));
d.loops(1).trafofn(2).args={}; 
d.loops(1).trafofn(1).fn=@(x,y) freqs(x(1))/baseratio;
d.loops(1).trafofn(1).args={};
% Turn off awg before changing clock rate.
clear pf;
pf(1).fn=@(x) smset('PulseLine',alloff);
pf(1).args={};
d.loops(1).prefn=[pf d.loops(1).prefn];

% Set AWG clock rate before running polarize function
clear pf;
pf(1).fn=@(x) smset(awgclk,freqs(1));
pf(1).args={};
d.loops(1).postfn = [ pf d.loops(1).postfn];

%%  
% Old file from Mikey do not run this.
%Change the mask to be smart about clock rate changes.
dbzmask=d.loops(1).prefn(3).args{2}{d.loops(1).prefn(3).args{1}(1)};
jmask=d.loops(1).prefn(3).args{2}{d.loops(1).prefn(3).args{1}(2)};
masks={dbzmask};
if 0  % Only fix T1
    for j=2:length(freqs)
        masks{j}=jmask*0;
    end
    for i=1:size(jmask,1)
        jfirst=find(jmask(i,:),1,'first');
        jlen=find(jmask(i,:),1,'last')-jfirst;
        for j=2:length(freqs)
            masks{j}(i,jfirst:(jfirst+ceil(jlen*freqs(j)/1e9)))=1;
        end
    end
else % Fix everything.
   ro=plsinfo('ro',d.data.pulsegroups(2).name);
   r=pdload('right');   
   for j=2:length(freqs)
       hwsr=freqs(j)/baseratio;
       masks{j}=false(size(jmask));
       nro=ro;
       nro(2)=(nro(2)-r.meas.time(2))*freqs(1)/freqs(j); % Actual start of readout.
       nro(2)=nro(2)+r.meas.time(2);       
       maskind = round(hwsr*cumsum(nro(2:3))*plsdata.tbase/freqs(1));
       %fprintf('Frequency %g: Readout %g-%g (%d-%d)/%d\n',freqs(j),nro(2:3),maskind(1:2),size(masks{j},2));
       masks{j}(1:2,maskind(1):maskind(2))=1;
   end
end
for j=2:length(freqs)
    masks{j}=logical(masks{j});
end
d.loops(1).prefn(3).args{1}=1:length(freqs);
d.loops(1).prefn(3).args{2}=masks;

% a bit of a hack because I'm not so clever
d.data.frange = drange;
d.data.freqs = freqs;

if 1 % more saving. scans are so long!
   d.saveloop = [2 1]; 
end

'Scan ready'