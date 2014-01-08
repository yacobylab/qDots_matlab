%% for crazy multi
rmpath('Z:\qDots\matlab\pulsecontrol2');
addpath('z:/qDots/matlab/pulsecontrol_sandbox/');
smdata.inst(14).cntrlfn = @smcATS660v4_v2;

%% regular

addpath('Z:\qDots\matlab\pulsecontrol2');
rmpath('z:/qDots/matlab/pulsecontrol_sandbox/');
smdata.inst(14).cntrlfn = @smcATS660v4;


%%
for ddt = [12]
bb=tic;
counter = 1;
offvalues =[0:50:300, 500:50:600, 800:50:900 1100:50:1200, 1500:50:1600];
%offvalues =[850:50:1050];% [0:50:200 500:50:600 1100:50:1200];
for offst = offvalues
nshts = 200; nrep = nshts; dt = ddt;
nevo = 50;
 
plslength=4; 

clear pg 
pg.pulses=12;
pg.ctrl='loop pack'; 
gen = '1';
pg.chan=[3 4]; pg.dict='right';

%params = [pulse length, max sep time, readdout time, sep time]
tev=offst+(0:(nevo-1))'; %50*ones(npls/2,1);
%tev = (1:2:2*nevo)';
%dt = 10;
%20ns steps
trep=(dt:dt:nrep*dt); trep = repmat(trep,1,nshts/length(trep))';
trep = flipud(trep);
%trep = [(10:10:250) 300:100:7700]';



pg.name=sprintf('dBz_rottest_2013_07_20_%s',upper(pg.dict(1)));pg.varpar=[trep; tev];%pg.varpar = repmat((dt:dt:50*dt),1,5)'; % 250, 4us This is the good one!
%pg.name=sprintf('dBz_rottest_2013_08_29_%s',upper(pg.dict(1))); pg.varpar=[trep; tev; tev; tev; tev];%pg.varpar = repmat((dt:dt:50*dt),1,5)';% 400, 6us
pg.params=[plslength max(pg.varpar) NaN 0]; 
%pg.name=sprintf('dBz_rottest_2013_07_23_%s',upper(pg.dict(1))); % 150, 8us
plsupdate(pg);
jgrpR=pg.name;

awgadd(pg.name);
awgcntrl('on start wait err raw');
%sm_setgradient; sm_getgradient
if sm_setgradient
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    d=smrun(fConfSeq2(pg.name,struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_%d_%d_R',offst,dt)));
    if any(isnan(d{1})); break; end
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end
fprintf('done with %d of %d in %d seconds \n\n',counter,length(offvalues),toc(bb));
counter = counter+1;

end
sleep

end
sleep

%% same scan as above, but forward time, and repeat estimation twice. 


bb=tic;
counter = 1;
offvalues =[0:50:300, 500:50:600, 850:50:900 1150:50:1200, 1550:50:1600 1950:50:2000];
dt=14; nrep=170; nevo = 50; plslength=4;
clear pg 
pg.pulses=12; pg.ctrl='loop pack'; gen = '1'; 
pg.chan=[3 4]; pg.dict='right';

trep=(dt:dt:nrep*dt)'; 
%trep=[trep; trep];
for offst = offvalues  
    %params = [pulse length, max sep time, readdout time, sep time]
    tev=offst+(0:(nevo-1))'; %50*ones(npls/2,1);        
    pg.varpar=[trep; tev; flipud(trep)];
    pg.params=[plslength max(pg.varpar) NaN 0];
    %dBz_rottest_2013_09_06 refers to 170 / 170 w/ 14 ns samp time. 
    %dBz_rottest_2013_09_07 refers to 170 / 130 w/ 14 ns samp time. 
    %dBz_rottest_2013_09_08 refers to 170 / 170 w/ 14 ns samp time and 0.9 us meas, 3.7 pulse. 
    %dBz_rottest_2013_09_08_2 refers to 170 / 170 w/ 14 ns samp time and 0.8 us meas, 300 ns delay, 3.6 pulse. 
    %dBz_rottest_2013_09_08_3 refers to 170 / 170 w/ 14 ns samp time and 0.7 us meas, 3.5 pulse. 
    %pg.name=sprintf('dBz_rottest_2013_09_08_3_%s',upper(pg.dict(1)));
    pg.name=sprintf('dBz_rottest_2013_09_08_4_%s',upper(pg.dict(1))); % 200-50-200, 4us
    try
        plsupdate(pg);
    catch 
        plsdefgrp(pg); 
    end
%     if offst==0
%         plsdefgrp(pg);
%     else
%         plsupdate(pg);
%     end
    jgrpR=pg.name;    
    awgadd(pg.name);
    awgcntrl('on start wait err raw');
    if sm_setgradient        
        d=smrun(fConfSeq2(pg.name,struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_%d_%d_R',offst,dt)));
        if any(isnan(d{1})); break; end                
    else
        fprintf('cannot lock gradient');
    end
    fprintf('done with %d of %d in %d seconds \n\n',counter,length(offvalues),toc(bb));
    counter = counter+1;    
end
sleep

%% Try out non uniform sampling


bb=tic;
counter = 1;
offvalues =[0:50:300, 500:50:600, 850:50:900 1150:50:1200, 1550:50:1600 1950:50:2000];
nevo = 50; 
clear pg 
pg.pulses=12; pg.ctrl='loop pack'; gen = '1'; 
pg.chan=[3 4]; pg.dict='right';
trand{1}=sort(floor(rand(100,1).*2400));
trand{2}=sort(floor(rand(200,1).*2400));
trand{3}=sort(floor(rand(50,1).*2400));
trand{5}=sort(floor(rand(100,1).*1900));
trand{4}=sort(floor(rand(50,1).*1400));
trand{6}=sort(floor(rand(100,1).*1400));
plslenmat=[4 4 4 3.5 3 3];
for j=1:length(trand) 
    for offst = offvalues
        %params = [pulse length, max sep time, readdout time, sep time]
        tev=offst+(0:(nevo-1))'; %50*ones(npls/2,1);
        pg.varpar=[trand{j}; tev];
        pg.params=[plslenmat(j) max(pg.varpar) NaN 0];
        pg.name=sprintf('dBz_rottest_2013_09_07_rand%d_%s_',j,upper(pg.dict(1)));
        try
            plsupdate(pg);
        catch
            plsdefgrp(pg);
        end
        jgrpR=pg.name;
        awgadd(pg.name);
        awgcntrl('on start wait err raw');
        if sm_setgradient
            d=smrun(fConfSeq2(pg.name,struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_%d_%d_R',offst,dt)));
            if any(isnan(d{1})); break; end
        else
            fprintf('cannot lock gradient');
        end
        fprintf('done with %d of %d in %d seconds \n\n',counter,length(offvalues),toc(bb));
        counter = counter+1;
    end
end
sleep
%% lets turn the scan on its side
nshts = 200; % number of points for estimation
rpeats = 1; %number of times to repepat each estimation point
nevo = 50; % number of evo times for measurement
evo = 1; % the evolution time
plslength=4; 
dt = 10; %time spacing for fourier estimate


clear pg 
pg.pulses=12;
pg.ctrl='loop pack'; 
gen = '1';
pg.chan=[3 4]; pg.dict='right';

for offst =0:28:2000
jgrpR={};
for evo = offst+(1:28)
tev=evo*ones(1,nevo)';
trep=(dt:dt:(nshts/rpeats)*dt); trep = repmat(trep,1,nshts/length(trep))';
pg.xval = evo;
%params = [pulse length, max sep time, readdout time, sep time]
pg.varpar=[trep; tev]; 
pg.params=[plslength max(pg.varpar) NaN 0]; %start sweep at eps=0. we dont want to avoid stp peak
%pg.nrep = ones(size(pg.varpar))'; pg.nrep(end) = nevo;
pg.name=sprintf('dBz_phase_2013_07_17c_%d_%s',length(jgrpR)+1,upper(pg.dict(1)));
%plsupdate(pg);
try
    plsupdate(pg)
catch
    plsdefgrp(pg);
end
jgrpR=[jgrpR, pg.name];
end
awgadd(jgrpR);
awgcntrl('on start wait err raw');
%sm_setgradient; sm_getgradient
if sm_setgradient
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    d=smrun(fConfSeq2([28:55],struct('nloop',1,'nrep',40, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_%d_R',offst)));
    if any(isnan(d{1})); break; end
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end

end
sleep



%% now lets try to estimate drift and difussion


for ddt = [12]
bb=tic;
counter = 1;
dvals = [0 100:100:1000]; % in microseconds
for dly = dvals;


nshts = 200; nrep = nshts; dt = ddt;

 
plslength=4;
clear pg 
pg.pulses=12;
pg.ctrl='loop pack'; 
pg.chan=[3 4]; pg.dict='right';

%params = [pulse length, max sep time, readdout time, sep time]


twait = dt*(1:floor(dly/plslength)); 
if max(twait) > 2400 %hack
   twait = [twait(1:100), twait(1:((length(twait))-100))]; 
end
trep=(dt:dt:nrep*dt); trep = repmat(trep,1,nshts/length(trep))';
pg.varpar=[trep; twait'; trep];


pg.name=sprintf('dBz_%d_2013_10_14_%s',length(pg.varpar),upper(pg.dict(1)));
pg.params=[plslength max(pg.varpar) NaN 0]; 

awgrm(27,'after'); awgclear('unused');
try 
    plsupdate(pg)
catch
    plsdefgrp(pg);
end
awgadd(pg.name);
awgcntrl('on start wait err raw');
%sm_setgradient; sm_getgradient
if sm_setgradient
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    d=smrun(fConfSeq2(pg.name,struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_%d_R',dly)));
    if any(isnan(d{1})); break; end
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end
fprintf('done with %d of %d in %d seconds \n\n',counter,length(dvals),toc(bb));
counter = counter+1;

end
sleep
end
sleep

%% group for just measuring over and over again


nshts = 200; nrep = nshts; dt = 12;
plslength=4;
clear pg 
pg.pulses=12;
pg.ctrl='loop pack'; 
pg.chan=[3 4]; pg.dict='right';

%params = [pulse length, max sep time, readdout time, sep time]


trep=(dt:dt:nrep*dt); trep = repmat(trep,1,nshts/length(trep))';
pg.varpar=[trep];


pg.name=sprintf('dBz_%d_2013_10_15_%s',length(pg.varpar),upper(pg.dict(1)));
pg.params=[plslength max(pg.varpar) NaN 0]; 

awgrm(27,'after'); awgclear('unused');
try 
    plsupdate(pg)
catch
    plsdefgrp(pg);
end
awgadd(pg.name);
awgcntrl('on start wait err raw');
gradtgt = abs(fbdata.gradtgt(2));
%sm_setgradient; sm_getgradient
if sm_setgradient
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    scan = fConfSeq2(pg.name,struct('nloop',256,'nrep',512, 'opts','raw pol logdbz','datachan','DAQ2'));
    scan.saveloop = [1 256];
    d=smrun(scan,smnext(sprintf('dBz_phase_raw_%.0fMHz_R',gradtgt)));
    if any(isnan(d{1})); break; end
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end
%fprintf('done with %d of %d in %d seconds \n\n',counter,length(dvals),toc(bb));
%counter = counter+1;

sleep


%% phase estimation with the extra markers for CDS

for ddt = [12]
bb=tic;
counter = 1;
offvalues =[0:50:300, 500:50:600, 800:50:900 1100:50:1200, 1500:50:1600];
%offvalues =[850:50:1050];% [0:50:200 500:50:600 1100:50:1200];
for offst = offvalues
nshts = 200; nrep = nshts; dt = ddt;
nevo = 50;
 
plen=4; 

clear pg 
pg.pulses=[71 12];
pg.ctrl='multi pack single'; 
gen = '1';
pg.chan=[3 4]; 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep'),'right'};

%params = [pulse length, max sep time, readdout time, sep time]
tev=offst+(0:(nevo-1))'; %50*ones(npls/2,1);
trep=((dt:dt:nrep*dt)-dt)';% trep = repmat(trep,1,nshts/length(trep))';
%trep = flipud(trep);

rtmp= pdload(pg.dict{end});
m_start = 1e-3*max([trep'])+sum(rtmp.reload.time)+.12;
m_dur = 1;
pg.name=sprintf('dBz_rot_CDS_2013_12_19_%s',upper(pg.dict{end}(1)));
pg.varpar={trep, tev};
pg.params = {[plen, max(trep),NaN,m_start,m_dur,0],[plen, max(tev), NaN,0]};

%pg.name=sprintf('dBz_rottest_2013_07_23_%s',upper(pg.dict(1))); % 150, 8us
plsupdate(pg);

awgadd(pg.name);
awgcntrl('on start wait err raw');
%sm_setgradient; sm_getgradient
if sm_setgradient
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    d=smrun(fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')),smnext(sprintf('dBz_phase_FPGA_%d_%d_R',offst,dt)));
    %d=smrun(fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',1, 'opts','pol logdbz','datachan','DAQ2')));
    if any(isnan(d{1})); break; end
    %fprintf('dBz = %.2f MHz\n',(sm_getgradient));
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end
fprintf('done with %d of %d in %d seconds \n\n',counter,length(offvalues),toc(bb));
counter = counter+1;

end
sleep

end
sleep


%% dBz phase estimation w/ the FPGA and vi control. 
dt = 12;
bb=tic;
counter = 1;
%offvalues =0:100:1700; % 2 ns sampling. 
%offvalues=1800:100:2400; 
offvalues=0:200:2000; %4 ns sampling
nreps=1024; nexp=length(offvalues); 
nsamps=199;
CS.MaxPulseCount.set(nsamps);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N elements',nsamps);
vi.Run([true]);
vi.SetControlValue('Reset',true);

for offst = offvalues
    nshts = nsamps+1; nevo = 50;    
    plen=4;
    clear pg; pg.pulses=[71 12];
    pg.ctrl='multi pack single';
    gen = '1'; pg.chan=[3 4];
    pg.dict={struct('prep','@adprep', 'read','@adread'),'right'};
    
    %params = [pulse length, max sep time, readdout time, sep time]
    tev=offst+4*(0:(nevo-1))'; %50*ones(npls/2,1);
    trep=((dt:dt:nshts*dt)-dt)';% trep = repmat(trep,1,nshts/length(trep))';
    %trep = flipud(trep);
    
    rtmp= pdload(pg.dict{end});
    m_start = 1e-3*max([trep'])+sum(rtmp.reload.time)+.12;
    m_dur = 1;
    pg.name=sprintf('dBz_rot_CDS_2013_12_19_%s',upper(pg.dict{end}(1)));
    pg.varpar={trep, tev};
    pg.params = {[plen, max(trep),NaN,m_start,m_dur,0],[plen, max(tev), NaN,0]};        
    plsupdate(pg);
    awgadd(pg.name);
    awgcntrl('on start wait err raw');    
    if sm_setgradient        
        scan=fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',nreps, 'opts','FPGA2 pol logdbz','datachan','DAQ2'));
        scan.data.FPGAThreshold=vi.GetControlValue('Threshold control'); 
        %d=smrun(scan,smnext(sprintf('dBz_phase_FPGA_%d_%d_R',offst,dt)));        
        d=smrun(scan); 
        if any(isnan(d{1})); break; end        
        fprintf(' \n \n');
    else
        fprintf('Cannot lock gradient');
    end
    fprintf('Done with %d of %d in %d seconds \n\n',counter,length(offvalues),toc(bb));
    counter = counter+1;
    
end
sleep

    
    