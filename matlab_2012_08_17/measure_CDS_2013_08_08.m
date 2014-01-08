%%
%Cells : (you can search these) 

% Initialize FPGA
% sweep a parameter
%FPGA gain sweep / offset sweep ramsey
%FPGA ramsey, time offset sweep
%FPGA Ramsey Ramsey 
%Rabi for pi/2 
%Make new CDS 
%ramsey-ramsey w/ multiple pulse times. 
%% Initialize FPGA
%Get application reference to Labview
e=actxserver('LabVIEW.Application');

%Open the FPGA project, and get reference to application instance of the
%project
projpath='Z:\John\Labview FPGA projects\delta Bz 1.6.2014\delta Bz project 1.6.2014.lvproj';
proj=invoke(e,'OpenProject',projpath);
f=proj.Application;

%Open the VI reference
vipath='Z:\John\Labview FPGA projects\delta Bz 1.6.2014\Read from CDS test data only 1.6.14.vi';
global vi;
vi=invoke(f,'GetVIReference',vipath);
vi.FPWinOpen=1;
vi.Run([true]);

%% offset sweep to find optimal working point on VCO
CS.DACGain.set(0);
ofvals = linspace(6,8,10);
d3=[];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
beta0=[0.016 0.002 150 0.01 500 0];
t2=zeros(1,length(ofvals)); 
scan=fConfSeq2_v2([28],struct('nloop',1,'nrep',512, 'opts','pol logdbz','datachan','DAQ2'));
for j = ofvals
   CS.DACOffset.set(j);
   if sm_setgradient
       sm_getgradient
       scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
       d=smrun(scan,smnext('Ramsey_CDS_off_R'));       
       sm_getgradient
       if any(isnan(d{1})); break; end;
       d3(end+1,:)=mean(d{1}(:,1001:end));
   else
       d3(end+1,:) = NaN(1,256);
      fprintf('no lock grad \n'); 
   end
   j
   smget('DMM')
   fprintf('done with %d of %d \n \n',size(d3,1),length(ofvals));
   figure(7); imagesc(4.*(1:500),ofvals(1:size(d3,1)),d3);
 
end
sleep
try
for j=1:length(ofvals)
    T=4.*(1:500);
    pars = fitwrap('plfit plinit',T,d3(j,:),beta0,fitfn, [1 1 1 1 1 0]);
    t2(j)=abs(pars(5));
    beta0=pars;
    figure(4);
    plot(ofvals(1:j),t2(1:j),'.');
end
end
sleep
%%
gainvals = round(linspace(1.3*2^22,2^23,32));
d3=[];
for j = gainvals
   CS.DACGain.set(j);
   if sm_setgradient
       sm_getgradient
       d=smrun(fConfSeq2_v2([31],struct('nloop',1,'nrep',256, 'opts','pol','datachan','DAQ2')));
       sm_getgradient
       if any(isnan(d{1})); break; end;
       d3(end+1,:)=mean(d{1}(:,1001:end));
   else
       d3(end+1,:) = NaN(1,500);
      fprintf('no lock grad \n'); 
   end
   figure(6); imagesc(1:500,gainvals(1:size(d3,1)),d3);
   fprintf('RV=%d,  DMM= %d \n',CS.ResultValue.get(),cell2mat(smget('DMM')))
   fprintf('G*RV + Offset = %d',(j/65536)*CS.ResultValue.get()+CS.DACOffset.get());
   fprintf('done with %d of %d \n \n',size(d3,1),length(gainvals));
end
sleep

%%
gainvals = round(linspace(1.5e6,2.5e6,32));
d3=[];
setpt= 6.7;
for j = gainvals
   CS.DACGain.set(j);
   sm_setgradient;
   for k = 1:30
       d=smrun(fConfSeq2_v2([31],struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2')));
       tmp2(k) = CS.ResultValue.get();
       tmp(k) = cell2mat(smget('DMM'));
   end
   ee = mean(tmp)-setpt;
   curr= CS.DACOffset.get();
   CS.DACOffset.set(curr-ee)
   if sm_setgradient
       %sm_getgradient
       scan = fConfSeq2_v2([31],struct('nloop',1,'nrep',256, 'opts','pol logdbz','datachan','DAQ2'));
       scan.data.CSOffset = CS.DACOffset.get();
       scan.data.CSGain = CS.DACGain.get();
       d=smrun(scan,smnext('ramsey_CDS_R'));
       %sm_getgradient
       if any(isnan(d{1})); break; end;
       d3(end+1,:)=mean(d{1}(:,1001:end));
   else
       d3(end+1,:) = NaN(1,500);
      fprintf('no lock grad \n'); 
   end
   figure(7); imagesc(1:500,gainvals(1:size(d3,1)),d3);
   fprintf('RV=%d,  DMM= %d \n',CS.ResultValue.get(),cell2mat(smget('DMM')))
   fprintf('done with %d of %d \n \n',size(d3,1),length(gainvals));
end
sleep


%%
%gainvals = round((4/3)*linspace(-2^23,2^23,32));


gainvals =linspace(0,-4e7,8);
%gainvals=0; 
setpt= 6.43;
scan2=(fConfSeq2_v2([28],struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2')));
scan2.loops(2).getchan = {'DMM'}; scan2.loops(2).npoints = 100;
d3=[];
scan=(fConfSeq2_v2([28],struct('nloop',1,'nrep',128, 'opts','pol logdbz','datachan','DAQ2')));
scan.loops.getchan={'DAQ2' 'Time' 'DMM'}; 
scan.loops.procfn(4)=scan.loops.procfn(3); 
scan.loops.procfn(3).dim=[]; scan.loops.procfn(3).fn=[];
scan.loops.procfn(1).fn(1).outchan=4; 
scan.disp(2).loop=1; scan.disp(2).channel=3; scan.disp(2).dim=1; 
scan.disp(1).loop=1; scan.disp(1).channel=1; scan.disp(1).dim=2; 
scan.loops.getchan{3}='DMM';
for j = gainvals
   CS.DACGain.set(j);
   sm_setgradient;
   d=smrun(scan2);
   ee =mean(d{end})-setpt;
   curr= CS.DACOffset.get();
   CS.DACOffset.set(curr-ee)
   if sm_setgradient
       %sm_getgradient
       %scan = fConfSeq2_v2([31],struct('nloop',1,'nrep',256, 'opts','pol logdbz','datachan','DAQ2'));
       scan.data.CSOffset = CS.DACOffset.get();
       scan.data.CSGain = CS.DACGain.get();
       scan.data.CSPreDelay=CS.PreDelay.get();
       scan.data.CSPostDelay=CS.PostDelay.get();
       scan.data.CSNSamps=CS.NSamps.get();
       scan.data.CSSampleCount=CS.SampleCount.get(); scan.data.setpt=setpt; 
       d=smrun(scan)%,smnext('ramsey_CDS_R'));
       %sm_getgradient
       if any(isnan(d{1})); break; end;
       d3(end+1,:)=mean(squeeze(d{1}(:,1001:end)));
   else
       d3(end+1,:) = NaN(1,500);
      fprintf('no lock grad \n'); 
   end
   figure(7); imagesc((0:499)*4,gainvals(1:size(d3,1)),d3);
   fprintf('RV=%d,  DMM= %d \n',CS.ResultValue.get(),cell2mat(smget('DMM')))
   fprintf('done with %d of %d \n \n',size(d3,1),length(gainvals));
end

fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
beta0=[0.016 0.002 75 0.01 500 0];
t2=zeros(1,length(gainvals));

for j=1:length(gainvals)
                 T=4.*(1:500);
                 pars = fitwrap('plfit plinit',T,d3(j,:),beta0,fitfn, [1 1 1 1 1 0]);
                 t2(j)=abs(pars(5)); 
                 beta0=pars; 
                 figure(4); 
                 plot(gainvals(1:j),t2(1:j),'.');    
end
sleep

%% Gain sweep
gainvals =[0 5.06e6];
setpt= CS.VCO_setpt; 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2([29],struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');
for j = gainvals
    CS.DACGain.set(j);
    CS.adjust_offset;
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        if exist('setpt','var')
            scan.data.setpt=setpt;
        end
        d=smrun(scan,smnext('ramsey_CDS_R'));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        d3(end+1,:) = NaN(1,500);  fprintf('Gradient would not lock \n');
    end
    figure(7); imagesc((0:499)*4,gainvals(1:size(d3,1)),d3);
    fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);
    fprintf('done with %i of %i scans \n\n',size(d3,1),length(gainvals));
end

%plot the t2*s, then the average dmm and daq val for each scan. these
%should be constant. 
beta0=[0.016 0.002 150 0.01 500 0];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
for j=1:length(gainvals)
    try
    T=4.*(0:499);
    pars = fitwrap('plfit plinit',T,d3(j,:),beta0,fitfn, [1 1 1 1 1 0]);        
    t2(j)=abs(pars(5));
    catch
       t2(j) = NaN; 
    end
    %beta0=pars;
    figure(4);
    plot(gainvals(1:j),t2(1:j),'.-');


end
figure(5); plot(gainvals,d2); 
figure(6); plot(gainvals,mean(dmm,2)); 
sleep

%% two scans
filenames{1}='ramsey_CDS_16db_R';
filenames{2}='ramsey_CDS_16db_0_R'; 
gainvals(1)= 5e6; 
gainvals(2)= 0; 

setpt= CS.VCO_setpt; 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2([28],struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');

%fix offset
for j=1:length(gainvals)
    CS.DACGain.set(gainvals(j)); 
    CS.adjust_offset; 
    
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        if exist('setpt','var')
            scan.data.setpt=setpt;
        end
        d=smrun(scan,smnext(filenames{j}));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        fprintf('Gradient would not lock \n');
    end
end
sleep
fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);


%% single scan

setpt= CS.VCO_setpt; 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2([29],struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');

%fix offset
CS.adjust_offset; 
    
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        if exist('setpt','var')
            scan.data.setpt=setpt;
        end
        d=smrun(scan,smnext('ramsey_CDS_R'));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        fprintf('Gradient would not lock \n');
    end
sleep
fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);






%% Sweep the MaxPulseCount
%Remember to set the gain and MaxPulseCount correctly first. 
%samps =linspace(500,990,75); samps=[samps 999];
CS.DACGain.set(6e6); 
CS.MaxPulseCount.set(850); 
%samps = [600 650 700 750 800 850 900 950 999];
samps=[1 50 100 150 200 250 300 350];
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2([28],struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');
orig_gain = CS.DACGain.get();
orig_samps = CS.MaxPulseCount.get();
for j = samps
    CS.MaxPulseCount.set(j);
    CS.DACGain.set(round(orig_gain*orig_samps/j));
    fprintf('setting gain to %d for MaxSampleCount = %d \n',CS.DACGain.get(),j);
    CS.adjust_offset()
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        scan.data.calibgain=CS.calibgain;
        scan.data.setpt=CS.VCO_setpt; scan.data.ampl=CS.ampl; scan.data.d_dbz=CS.d_dbz;         
        d=smrun(scan,smnext(sprintf('ramsey_CDS_%d_R',j)));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        d3(end+1,:) = NaN(1,500);  fprintf('Gradient would not lock \n');
    end
    figure(7); imagesc((0:499)*4,samps(1:size(d3,1)),d3);
    fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);
end
sleep
%% measure dbz with the CDS
CS.DACGain.set(4.3e4); 
sm_setgradient; sm_getgradient
DAQvals = [];
CDSvals = [];
for j = 1:128
    d=smrun(fConfSeq2_v2([29],struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2')));
    DAQvals(j,:) = squeeze(d{1});
    CDSvals(j,:) = CS.getBuffers(100); % takes about 4s
    fprintf('done with %d \n',j);
end
sleep
cv = CDSvals(:,1:99); dv = DAQvals(:,2:100);
figure(2)
plot(dv(:),cv(:),'.')


%% check buffer vals to compare to DMM/DAQ 

sm_setgradient; sm_getgradient
  DAQvals = [];
  CDSvals = [];
  DMMvals=[];        
  conf.offset=CS.DACOffset.get(); conf.gain=CS.DACGain.get(); conf.pre=CS.PreDelay.get(); 
  conf.post=CS.PostDelay.get(); conf.sample=CS.SampleCount.get(); conf.nsamps=CS.NSamps.get();
        for j = 1:128
            d=smrun(fConfSeq2_v2([28],struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2')));
            DAQvals(j,:) = squeeze(d{1});
            CDSvals(j,:) = CS.getBuffers(999); % takes about 4s
            DMMvals(j)=cell2mat(smget('DMM'));
            fprintf('done with %d \n',j);
        end
        sleep
        %filename=sprintf('BuffDAQ_%02d',17); 
        %save(filename,'DAQvals','CDSvals','DMMvals','conf')

%% check measurement via DMm v. DAQ. sweep post/predelay
scan=fConfSeq2_v2([28],struct('nloop',1,'nrep',512, 'opts','pol','datachan','DAQ2')); 
scan.loops.getchan={'DAQ2' 'Time' 'DMM'}; 
scan.loops.procfn(4)=scan.loops.procfn(3); 
scan.loops.procfn(3).dim=[]; scan.loops.procfn(3).fn=[];
scan.loops.procfn(1).fn(1).outchan=4; 
scan.disp(2).loop=1; scan.disp(2).channel=3; scan.disp(2).dim=1; 
scan.disp(1).loop=1; scan.disp(1).channel=1; scan.disp(1).dim=2; 
clear cv; clear dv; 
figure(2); clf;

for i=1:4
    CS.PostDelay.set(i+1);
    for k=1:4;     
        CS.PreDelay.set(k+2); 
        sm_setgradient; sm_getgradient   
        d=smrun(scan);
        DAQvals = squeeze(d{1});
        DMMvals=squeeze(d{3}); 
        cv{i,k}=DMMvals; dv{i,k} = mean(DAQvals(:,2:1000),2);
        figure(2);
        subplot(4,4,sub2ind([4 4],k,i)); 
        plot(dv{i,k}(:),cv{i,k}(:),'.'); 
    end

end
       filename=sprintf('DMMDAQ_%02d',2); 
        save(filename,'cv','dv')
        sleep
        
%% check measurement via DMM v. DAQ. Keep max time const. Sweep Sc/PostDel
scan=fConfSeq2_v2([28],struct('nloop',1,'nrep',512, 'opts','pol','datachan','DAQ2')); 
scan=addchan(scan,'DMM');
clear cv; clear dv; 
figure(2); clf;
gain=CS.DACGain.set(6e6); 
sc=CS.SampleCount.set(2); 
CS.adjust_offset; 

for i=1:7
    CS.PostDelay.set(i);    
    CS.SampleCount.set(8-i);
    CS.DACGain.set(gain.*2./(8-i)); 
    if sm_setgradient
        sm_getgradient   
        d=smrun(scan);
        DAQvals = squeeze(d{1});
        DMMvals=squeeze(d{3}); 
        cv{i}=DMMvals; dv{i} = mean(DAQvals(:,2:1000),2);
        figure(2);
        subplot(1,7,i)
        plot(dv{i}(:),cv{i}(:),'.');          
    end

end
       filename=sprintf('DMMDAQ_%02d',6); 
        save(filename,'cv','dv')
        sleep
        
        
%% run normal dBZ w/ the box 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2([28],struct('nloop',1,'nrep',512, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');

if sm_setgradient
    scan.data.CSOffset = CS.DACOffset.get();
    scan.data.CSGain = CS.DACGain.get();
    scan.data.CSPreDelay=CS.PreDelay.get();
    scan.data.CSPostDelay=CS.PostDelay.get();
    scan.data.CSNSamps=CS.NSamps.get();
    scan.data.CSSampleCount=CS.SampleCount.get();
    d=smrun(scan,smnext('dBz_CDS_R'));
    if any(isnan(d{1})); break; end;
    d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
    d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
    dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
else
    fprintf('Gradient would not lock \n');
end

sleep

%% new type of offset scan, to compare to a VCO sweep. 

offvals=linspace(5,8,32); 
scan = fConfSeq2_v2([31],struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2'));
scan2 = fConfSeq2_v2([30],struct('nloop',100,'nrep',50, 'opts','logdbz pol offok','datachan','DAQ2'));
CS.DACGain.set(0); 

for j=offvals
    CS.DACOffset.set(j);
    %sm_setgradient; sm_getgradient;
    smrun(scan);
    smget('DMM')
    if sm_setgradient; %sm_getgradient; 
        d=smrun(scan2,smnext(sprintf('rotramsey_VCO_%.2f_R',j)));
        if any(isnan(d{1})); break; break; end;
    end
end

sleep

%% new type of offset scan, to compare to a VCO sweep (uses real ramsey in rotating frame)

offvals=linspace(6,8,16); 
scan = fConfSeq2_v2('rotramsey_CDS_2013_07_28_R2',struct('nloop',1,'nrep',1, 'opts','pol','datachan','DAQ2'));
scan2 = fConfSeq2_v2('dBz_rot_ramsey_2013_08_27',struct('nloop',50,'nrep',200, 'opts','logdbz pol offok','datachan','DAQ2'));
CS.DACGain.set(0); 

for j=offvals
    CS.DACOffset.set(j);
    %sm_setgradient; sm_getgradient;
    smrun(scan);
    smget('DMM')
    if sm_setgradient; %sm_getgradient; 
        d=smrun(scan2,smnext(sprintf('dbz_rotramsey_VCO_%.2f_R',j)));
        if any(isnan(d{1})); break; break; end;
    end
end

sleep
%% new gain scan to go with above
gainvals = [0 linspace(-1e6,-4e6,20)];%linspace(1e6, 1e7,32)];
setpt= CS.VCO_setpt; 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2('dBz_rot_ramsey_CDS_2013_08_28',struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');
for j = gainvals
    CS.DACGain.set(j);
    fprintf('setting gain to %02d \n',j);
    CS.adjust_offset;
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        if exist('setpt','var')
            scan.data.setpt=setpt;
        end
        d=smrun(scan,smnext('ramsey_CDS_R'));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        d3(end+1,:) = NaN(1,400);  fprintf('Gradient would not lock \n');
    end
    CS.NSamps.get()
    figure(7); imagesc((0:399)*2,gainvals(1:size(d3,1)),d3);
    fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);
    fprintf('done with %i of %i scans \n\n',size(d3,1),length(gainvals));
end

figure(5); plot(gainvals,d2); 
figure(6); plot(gainvals,mean(dmm,2)); 
sleep

%% gain scan, 

gainvals = [linspace(-1e7,-9e7,8)]%[0 linspace(5e6,5e7,8)];%linspace(1e6, 1e7,32)];
setpt= CS.VCO_setpt; 
d3=[]; dmm=[]; d2=[]; 
scan=(fConfSeq2_v2('dBz_rot_ramsey_CDS_2013_08_28',struct('nloop',1,'nrep',1024, 'opts','pol logdbz','datachan','DAQ2')));
scan=addchan(scan,'DMM');
for j = gainvals
    %fprintf('setting gain to %f \n',j);
    CS.DACGain.set(j);
    CS.adjust_offset;
    if sm_setgradient
        scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
        if exist('setpt','var')
            scan.data.setpt=setpt;
        end
        d=smrun(scan,smnext('ramsey_CDS_R'));
        if any(isnan(d{1})); break; end;
        d3(end+1,:)=mean(squeeze(d{1}(:,1001:end))); %d3 stores oscillations
        d2(end+1)=mean(mean(squeeze(d{1}(1:1000)))); %d2 stores shots
        dmm(end+1,:)=squeeze(d{3}); %dmm has the dmm voltage for each line.
    else
        d3(end+1,:) = NaN(1,400);  fprintf('Gradient would not lock \n');
    end
    CS.NSamps.get()
    figure(7); imagesc((0:399)*2,gainvals(1:size(d3,1)),d3);
    fprintf('Average DMM: %1.3f V. SetPt: %1.3f \n',mean(d{3}),setpt);
    fprintf('done with %i of %i scans \n\n',size(d3,1),length(gainvals));
end

figure(5); plot(gainvals,d2); 
figure(6); plot(gainvals,mean(dmm,2)); 
sleep
%% sweep a parameter
parname = 'SampleCount';
parvals = 5:12;

results = struct();
for j = 1:length(parvals)
    CS.(parname).set(parvals(j));
    d=smrun(fConfSeq2_v2('dBz_rot_ramsey_CDS_2013_08_28',struct('nloop',1, 'nrep',1,'opts','','datachan','DAQ2')));
    if CS.MaxPulseCount.get() ~=CS.NSamps.get()
        fprintf('missed triggers\n')
        keyboard
        break
    else
    results(j).DAQ = d{1}(1:1000);
    results(j).CDS = CS.getBuffers(CS.NSamps.get());
    results(j).(parname) = parvals(j);
    end
end
sleep
n_horiz = ceil(length(parvals)/2);
n_vert = 2;
figure(4); clf; 
for j = 1:length(results)
    subplot(n_vert,n_horiz,j);
    plot(results(j).DAQ(2:991),results(j).CDS,'.');
    tmp = corrcoef(results(j).DAQ(2:991),results(j).CDS);
    title(sprintf('%s = %d\n r^2=%.2f',parname,results(j).(parname),abs(tmp(2,1))))
    xlabel('DAQ');
    ylabel('CDS');
end

figure(5); clf; 
for j = 1:length(results)
    subplot(n_vert,n_horiz,j);
    [h, b] = hist(results(j).CDS,50);
    %[h2, b2] = hist(results(j).DAQ(2:991),50);
    %plotyy(b,h,b2,h2);
    plot(b,h);
    title(sprintf('%s = %d\n r^2=%.2f',parname,results(j).(parname),abs(tmp(2,1))))
    xlabel('CDS Readout');
    ylabel('Counts');
end

%% grab fpga data and add to scan

nrep = 512*2;
cc = struct('file_num',0,'opts','freqs array'); %can add 'file_num',sum_num
%plsgrp = 'dBz_rot_CDS_2013_12_19_R';
plsgrp = 'dBz_rot_ramsey_CDS_2013_12_18';
scan = fConfSeq2_v2(plsgrp,struct('nloop',1,'nrep',nrep, 'opts','pol logdbz','datachan','DAQ2'));
scan.cleanupfn(end+1).fn = @get_FPGA_freqs;
scan.cleanupfn(end).args = {cc};


%% single sweep of times used to verify
%FPGA, time offset sweep 
offvals =[50, 500:50:600, 800:50:900 1100:50:1200, 1500:50:160];
dt = 12; plen = 4;
bb=tic; counter = 1;
mean_out0=17975
mean_freq=144; 
vi.SetControlValue('Gain',0); 
vi.SetControlValue('Offset',mean_out0); 
nreps=1024; nexp=length(offvalues);  
nsamps=199; 
CS.MaxPulseCount.set(nsamps); 
    vi.SetControlValue('N exps',nexp); 
    vi.SetControlValue('N reps',nreps); 
    vi.SetControlValue('N elements',nsamps);    
    try
    vi.Run([true]);
    catch
    end


for offst = offvals
clear pg 
pg.pulses=[71 72];
pg.ctrl='multi pack single'; 
pg.chan = [3 4];
rtmp = pdload('right');
vp1 = (dt:dt:200*dt)'-dt;
vp2 =offst+(0:49)';
m_start = 1e-3*max([vp1'])+sum(rtmp.reload.time)+.12;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
pars1 = [plen, max(vp1), NaN, m_start, m_dur, 0]; pars2=[plen, 5.8, 0];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep'),'right'};
pg.varpar={vp1, vp2};%128+0*(1:128)';
pg.params = {pars1, pars2};
pg.name=sprintf(('dBz_rot_ramsey_CDS_2013_12_27'),upper(pg.dict{end}(1)));
plsupdate(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');

%now make the scan
scan = fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',nreps, 'opts','pol logdbz FPGA2','datachan','DAQ2'));

if sm_setgradient    
    d=smrun(scan,smnext(sprintf('rot_ramsey_FPGA_%d_%d_R',offst,dt)));    
    if any(isnan(d{1})); break; end    
    fprintf(' \n \n');
else
   fprintf('cannot lock gradient'); 
end
fprintf('done with %d of %d in %d seconds \n\n',counter,length(offvals),toc(bb));
counter = counter+1;

end
 vi.SetControlValue('stop',true)
sleep

%% Run a gain sweep using the FPGA 
%FPGA gain sweep / offset sweep ramsey
%to change this to an offset scan, sweep mean_out and set offset
%accordingly. 
%info for dummies: VCO is about 6 Mhz / V, so 1 Mhz = 0.16 V= 500 bits. 
% mean_freq=144; offset0=7895; g0=70.05; 
% mean_freq2=350;
mean_out0=18350; mean_freq=67; 

%mean_out0=offset0+g0*mean_freq; 
nreps=1024; nexp=9; nsamps=199; 
CS.MaxPulseCount.set(nsamps); 
    vi.SetControlValue('N exps',nexp); 
    vi.SetControlValue('N reps',nreps); 
    vi.SetControlValue('N elements',nsamps);    
    vi.Run([true]);
    vi.SetControlValue('Reset',true);

if 1    
clear pg 
pg.pulses=[71 72]; dt=12; 
pg.ctrl='multi pack single'; 
pg.chan = [3 4]; plen=4; 
rtmp = pdload('right');
vp1 = (dt:dt:200*dt)'-dt;
vp2 =16*(0:65)';
m_start = 1e-3*max([vp1'])+sum(rtmp.reload.time)+.12;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
pars1 = [plen, max(vp1), NaN, m_start, m_dur, 0]; pars2=[plen, norm(rtmp.sep.val)/sqrt(2), 0];
etmp=rtmp.exch; etmp.val(1:3) =[rtmp.markerburst.val(1:2),norm(rtmp.sep.val)/sqrt(2)]; %exchange at sep. 
% adreadtmp=rtmp.adread; adreadtmp.val(3)=-8/6;
% adpreptmp=rtmp.adprep(1); adpreptmp.val(3)=-8/6;
% RFypreptmp=rtmp.RFyprep;
% RFypreptmp.val(1)=-8/6; RFypreptmp.val(3)=6;
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
%pg.dict={struct('prep',adpreptmp, 'read2',adreadtmp,'prep2',RFypreptmp,'read',RFypreptmp,'exch',etmp),'right'};
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep','exch',etmp),'right'};
pg.varpar={vp1, vp2};%128+0*(1:128)';
pg.params = {pars1, pars2};
pg.name=sprintf(('dBz_rot_ramsey_CDS_2013_01_01'),upper(pg.dict{end}(1)));
plsupdate(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');
end
scan=fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA2 logdbz pol','datachan','DAQ2','FPGA','freqs array'));    
%    mean_outvals=linspace(mean_out0-4000,mean_out0+6000,nexp)
gainvals=[0 linspace(36,39,nexp-1)]; 
%gainvals=[37.25];
%gain=35; 
 for j=1:length(gainvals)
     gain=gainvals(j); 

% for j=1:length(mean_outvals)
%     mean_out=mean_outvals(j)
    offst=mean_out0-gain*mean_freq2;     
    %offst=mean_out-gain*mean_freq2;     
    vi.SetControlValue('Offset',offst); 
    vi.SetControlValue('Gain',gain);     
    if sm_setgradient
        smrun(scan,smnext('rot_ramsey_FPGA_R'));
    end
end
 vi.SetControlValue('stop',true)
sleep

%% Rabi Sweep to find pi/2 time. 
%Rabi for pi/2 
if 1 
clear pg 
pg.pulses=[91]; pg.ctrl='pack'; 
pg.chan = [3 4]; plen=4; 
rtmp = pdload('right');
vp2 =2*(0:65)';
pars2=[plen 0];
% adreadtmp=rtmp.adread; adreadtmp.val(3)=-8/6; 
% adpreptmp=rtmp.adprep(1); adpreptmp.val(3)=-8/6; 
 RFypreptmp=rtmp.RFyprep; 
 RFypreptmp.val(1)=-1; RFypreptmp.val(3)=5.8; 
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep','read','@adread','prep2',RFypreptmp),'right'};
pg.varpar=vp2;%128+0*(1:128)';
pg.params = pars2;
pg.name=sprintf(('rabi_2014_01_02'),upper(pg.dict{end}(1)));
plsupdate(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');

end

mean_out0=18350;
vi.SetControlValue('Offset',mean_out0);
vi.SetControlValue('Gain',0);
try
    vi.Run([true]);
catch
end
scan=fConfSeq2_v2(pg.name,struct('nloop',200,'nrep',64,'opts','logdbz pol','datachan','DAQ2'));
scan.loops(1).setchan='RFfreq3'; 
scan.loops(1).rng=[160e6 220e6];
scan.loops(2).npoints=10;
scan.loops(2).setchan='count';
scan.loops(2).rng = [1, scan.loops(2).npoints];
if sm_setgradient 
    smrun(scan,smnext('rabi_R'));
end


%% ramseyramsey!
%FPGA Ramsey Ramsey 
%estimate detuning w/ ramsey then do ramsey. 
%currently uses multiple time offsets. 
if 1

clear pg 
dt=16; num_pulse=4; num_evo=floor(120/num_pulse); %assumes we have 120 total evos. 
%num_evo=30;
plen_evo=5;
offvalstm=dt*num_evo*[0:num_pulse-1];
pg.pulses=[90 72]; 
pg.ctrl='multi single'; 
pg.chan = [3 4]; plen=4; 
rtmp = pdload('right'); 
%nsamps=60;
nsamps=50; 
dt2=16; 
vp1 = (dt2:dt2:(nsamps+1)*dt2)'-dt2;

%vp1=[vp1; vp1(2:end)]; %used for double sampling. 
nests=length(vp1)-1;
%vp1=sort(vp1); %used to go in order. 
%nsamps=length(vp1)-1;
%vp2=flipud(vp2); 
m_start = sum(rtmp.reload.time)+.12+sum([rtmp.adprep.time])+rtmp.adread.time;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
pars1 = [plen, norm(rtmp.sep.val)/sqrt(2), m_start, m_dur, 0]; pars2=[plen_evo, norm(rtmp.sep.val)/sqrt(2), 0];
etmp=rtmp.exch; etmp.val(1:3) =[rtmp.markerburst.val(1:2),norm(rtmp.sep.val)/sqrt(2)];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep','exch',etmp),'right'};
pg.params = {pars1, pars2};
%namepat='ramsey_rot_ramsey_CDS_2014_01_06_6_%s_%d'; 
namepat='ramsey_rot_ramsey_CDS_2014_01_08_2_%s_%d'; 
burnpts=zeros(5,1); 
for k=1:length(offvalstm)
    vp2 =offvalstm(k)+dt*[burnpts; (0:num_evo-1)']; %burn the first 5 points to let the FPGA update.
    pg.varpar={vp1, vp2};%128+0*(1:128)';
    pg.name=sprintf(namepat,upper(pg.dict{end}(1)),offvalstm(k));

    try
        plsupdate(pg);
    catch
        plsdefgrp(pg);
    end

    awgadd(pg.name);
end
awgcntrl('on start wait err raw');
end

mean_out0=17500; mean_freq=80; 
%cntrl_out=18350;
cntrl_out=17500; 
nreps=1024; nexp=4;
CS.MaxPulseCount.set(nests);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N index',nsamps);
vi.SetControlValue('N est',nests);
vi.SetControlValue('control output',cntrl_out);
vi.SetControlValue('phase',pi); %phase shift of pi for ramsey
%vi.SetControlValue('phase',2.6426); 

try
    vi.Run([true]);
catch
end

vi.SetControlValue('Reset',true);

for k=1:num_pulse
%for k=1
    
%scan=fConfSeq2_v4(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));    
pls_name=sprintf(namepat,upper(pg.dict{end}(1)),offvalstm(k)); %check that this is set correctly 
scan=fConfSeq2_v4(pls_name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));    
%gainvals=[linspace(27,29,nexp-1)];
gainvals=0;
%outvals=linspace(20250,20750,nexp-1);
for j=1:length(gainvals)
%for j=1:length(outvals)
    %gain=gainvals(1); 
    %offst=outvals(j)-gain*mean_freq;     
    gain=gainvals(j);
    offst=mean_out0-gain*mean_freq
    vi.SetControlValue('Offset',offst); 
    vi.SetControlValue('Gain',gain);     
    vi.SetControlValue('estimation dt',dt2);   
    
    if sm_setgradient
        %d=smrun(scan,smnext('ramsey_rot_ramsey_FPGA_R'));
        d=smrun(scan,smnext(sprintf('ramsey_rot_ramsey_FPGA_R_%d_%d',j,offvalstm(k))));
        %d=smrun(scan)
    end
    if any(isnan(d{1}))
       break 
    end
end
    if any(isnan(d{1}))
       break 
    end

end
 vi.SetControlValue('stop',true)
sleep

%% Sweep the number of estimations to see what the optimum is
%FPGA Ramsey Ramsey 
%estimate detuning w/ ramsey then do ramsey. 

ests=95%[75 100 125 150 175 200 225 250]
ngp = length(awgdata.pulsegroups);
brk = 0;
for j=1:length(ests)

nests=ests(j);
nsamps=nests;
clear pg 
est_dt=16; 
evo_dt=1;
num_evo=130; 
pg.pulses=[90 72]; 
pg.ctrl='multi pack single'; 
pg.chan = [3 4]; 
est_plen=5;
evo_plen=5;
rtmp = pdload('right'); 
vp1 = (est_dt:est_dt:(nests+1)*est_dt)'-est_dt;
vp2=(evo_dt:evo_dt:(num_evo)*evo_dt)-evo_dt;
vp2=[0 0 0 0 0 vp2]';
m_start = sum(rtmp.reload.time)+.12+sum([rtmp.adprep.time])+rtmp.adread.time;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
pars1 = [est_plen, norm(rtmp.sep.val)/sqrt(2), m_start, m_dur, 0]; pars2=[evo_plen, norm(rtmp.sep.val)/sqrt(2), 0];
etmp=rtmp.exch; etmp.val(1:3) =[rtmp.markerburst.val(1:2),norm(rtmp.sep.val)/sqrt(2)];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep','exch',etmp),'right'};
pg.params = {pars1, pars2};
namepat='ramsey_rot_ramsey_CDS_2014_01_06_%d'; 

pg.varpar={vp1, vp2};%128+0*(1:128)';
pg.name=sprintf(namepat,nests);   

try
plsupdate(pg);
catch
plsdefgrp(pg);
end

awgrm(ngp,'after');
awgclear('unused'); awgclear('unused');
awgadd(pg.name);    

awgcntrl('on start wait err raw');

mean_out0=18350; mean_freq=67; 
nreps=1024; nexp=5;
CS.MaxPulseCount.set(nests);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N index',nsamps);
vi.SetControlValue('N est',nests);
vi.SetControlValue('control output',mean_out0);

try
    vi.Run([true]);
catch
end

vi.SetControlValue('Reset',true);

scan=fConfSeq2_v4(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));    
scan.data.fourier_inds=[2:nests+1];
scan.data.evo_inds=[nests+2:nests+num_evo+1]+5;
%gainvals=[0 linspace(55,57,nexp-1)];
gainvals=56; 
for jj=1:length(gainvals)
    gain=gainvals(jj); 
    offst=mean_out0-gain*mean_freq;     
    vi.SetControlValue('Offset',offst); 
    vi.SetControlValue('Gain',gain);     
    
    if sm_setgradient
        %d=smrun(scan,smnext('ramsey_rot_ramsey_FPGA_R'));
        d=smrun(scan,smnext(sprintf('ramsey_rot_ramsey_FPGA_R_%d',nests)));
        %d=smrun(scan)
    end
    if any(isnan(d{1}))
       brk = 1;
        break; break;
    end
end

if brk; break; end
end

 vi.SetControlValue('stop',true)
sleep

%% ramsey-ramsey w/ multiple pulse times. 
%let's try pulses lengths 4, 5, 6. start by doing it in one pulse. 
% appears to be 2.89 us of overhead here. so (0:69), (70:130), (131:194) 
%start w/ everything hardcoded. 
mean_out0=18350; mean_freq=67; 
nreps=1024; nexp=6; nsamps=130; nests=nsamps; 
CS.MaxPulseCount.set(nests);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N index',nsamps);
vi.SetControlValue('N est',nests);
vi.SetControlValue('control output',mean_out0);

try
    vi.Run([true]);
catch
end

vi.SetControlValue('Reset',true);


clear pg 
dt=16; %num_evo=43; num_pulse=floor(196/num_evo);
%offvalstm=dt*num_evo*[0:num_pulse-1];
pg.pulses=[90,90,72]; 
pg.ctrl='multi single'; 
pg.chan = [3 4]; 
rtmp = pdload('right'); 
nsamps=130; 
vp1 = {dt*[0:68]',dt*[69:130]'};
nests=nsamps;

m_start = sum(rtmp.reload.time)+.12+sum([rtmp.adprep.time])+rtmp.adread.time;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pars1 = [0, norm(rtmp.sep.val)/sqrt(2), m_start, m_dur, 0]; pars2=[0, norm(rtmp.sep.val)/sqrt(2), 0]; %plength is dummy here, set later in cell
pars_est={pars1, pars1}; pars_est{1}(1)=4; pars_est{2}(1)=5;
pars_evo={pars2}; 
etmp=rtmp.exch; etmp.val(1:3) =[rtmp.markerburst.val(1:2),norm(rtmp.sep.val)/sqrt(2)];

pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep','exch',etmp),'right'};
namepat='ramsey_rot_ramsey_CDS_2014_01_06_2_%s';
%vp2set= {dt*[0:68]',dt*[69:130]',dt*[131:193]'};
vps=dt*[0:193]';
% max evo time for 4 us: 68*16, 5 us: 130 * 16 6 us: 193 * 16 
vp2set={vps(1:32),vps(33:64),vps(65:96),vps(97:128),vps(129:160),vps(161:192)};
plenset=[4;4;5;5;6;6];
burnpts=zeros(5,1);


ngp = length(awgdata.pulsegroups);
%ngp=29
%pls_name=sprintf(namepat,upper(pg.dict{end}(1)),offvalstm(k)); %check that this is set correctly 
%scan=fConfSeq2_v4(pls_name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));    
%gainvals=[0 linspace(55,57,nexp-1)];
gainvals=56; 
%for k=1:num_pulse

for k=1:1%length(vp2set)
    vp2 = vp2set{k}; 
    pg.name=sprintf(namepat,upper(pg.dict{end}(1)));
    vp2=[burnpts; vp2];
    pars_evo{1}(1)=plenset(k); 
    pg.params = {pars_est{:}, pars_evo{:}};
    pg.varpar={vp1{:}, vp2};%128+0*(1:128)';    
    awgrm(ngp,'after');
    awgclear('unused');
    try 
        plsupdate(pg);
    catch
        plsdefgrp(pg);
    end
    awgadd(pg.name);
    awgcntrl('on start wait err raw');
    %scan=fConfSeq2_v4(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));    
    scan=fConfSeq2_v4(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA logdbz pol','datachan','DAQ2'));
    %gainvals=[0 linspace(55,57,nexp-1)];
  for j=1:length(gainvals)
        gain=gainvals(j);
        offst=mean_out0-gain*mean_freq;
        vi.SetControlValue('Offset',offst);
        vi.SetControlValue('Gain',gain);
        
        if sm_setgradient
            %d=smrun(scan,smnext('ramsey_rot_ramsey_FPGA_R'));
            d=smrun(scan,smnext('ramsey_rot_ramsey_FPGA_R_multtms'));
            %d=smrun(scan)
        end
        if any(isnan(d{1}))
            break
        end
    end
end
 vi.SetControlValue('stop',true)
sleep



%% Single scan 
mean_freq=144; offset0=7895; 
g0=70.05; 
mean_out=offset0+g0*mean_freq; 
nreps=1024; nexp=1;
nsamps=199; 
vi.SetControlValue('Gain',0);
vi.SetControlValue('Offset',mean_out);
CS.MaxPulseCount.set(nsamps);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N elements',nsamps);
vi.Run([true]);

clear pg 
pg.pulses=[71 72];
pg.ctrl='multi pack single'; 
pg.chan = [3 4];
rtmp = pdload('right');
vp1 = (dt:dt:200*dt)'-dt;
vp2 =20*(0:49)';
m_start = 1e-3*max([vp1'])+sum(rtmp.reload.time)+.12;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
pars1 = [plen, max(vp1), NaN, m_start, m_dur, 0]; pars2=[plen, 5.8, 0];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep'),'right'};
pg.varpar={vp1, vp2};%128+0*(1:128)';
pg.params = {pars1, pars2};
pg.name=sprintf(('dBz_rot_ramsey_CDS_2013_12_27'),upper(pg.dict{end}(1)));
plsupdate(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');

scan=fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA2 logdbz pol','datachan','DAQ2'));
if sm_setgradient
    smrun(scan,smnext('rot_ramsey_FPGA_R'));
end

 vi.SetControlValue('stop',true)
sleep

%% Single scan 
%defunct. ramsey ramsey w/ just one pulse 
mean_freq=144; offset0=7895; 
g0=70.05; 
mean_out=offset0+g0*mean_freq; 
nreps=1024; nexp=1;
nsamps=49; 
vi.SetControlValue('Gain',0);
vi.SetControlValue('Offset',mean_out);
CS.MaxPulseCount.set(nsamps);
vi.SetControlValue('N exps',nexp);
vi.SetControlValue('N reps',nreps);
vi.SetControlValue('N elements',nsamps);
vi.Run([true]);

dt=20; 
clear pg 
pg.pulses=[72]; 
pg.ctrl='pack single'; 
pg.chan = [3 4];
rtmp = pdload('right'); 
vp1 = (dt:dt:50*dt)'-dt;
vp2 =20*(0:49)';
m_start = 1e-3*max([vp1'])+sum(rtmp.reload.time)+.12;
m_dur = 1;%rtmp.meas.time(1:3)*[1;-1;-1];
%pars1 = [plen, max(vp1), NaN, m_start, m_dur, 0]; 
pars2=[plen, 5.8, 0];
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pg.dict={struct('prep','@adprep', 'read2','@adread','prep2','@RFyprep','read','@RFyprep'),'right'};
pg.varpar=[vp1; vp2];%128+0*(1:128)';
assert(nsamps ==length(vp1)-1);
pg.params = pars2;
pg.name=sprintf(('ramsey_rot_ramsey_CDS_2014_01_01'),upper(pg.dict{end}(1)));
plsdefgrp(pg);
awgadd(pg.name);
awgcntrl('on start wait err raw');

scan=fConfSeq2_v2(pg.name,struct('nloop',1,'nrep',nreps,'opts','FPGA2 logdbz pol','datachan','DAQ2'));
if sm_setgradient
    smrun(scan,smnext('rot_ramsey_FPGA_R'));
end

 vi.SetControlValue('stop',true)
sleep

%% offset scan
trick_scan = fConfSeq2_v2('dBz_rot_ramsey_CDS_2013_08_28',struct('nloop',1, 'nrep',1,'opts','','datachan','DAQ2'));
real_scan = fConfSeq2_v2('dBz_rot_ramsey_2013_08_27',struct('nloop',100, 'nrep',50,'opts','pol logdbz','datachan','DAQ2'));
offvals = linspace(5,8,32);
npls = 201;

d3 = zeros(length(offvals),npls);
CS.DACGain.set(0);
for j = 1:length(offvals)
    fprintf('setting offset to %.2d \n',offvals(j));
    CS.DACOffset.set(offvals(j));
    smrun(trick_scan);
    if sm_setgradient
       sm_getgradient
       scan.data.CSOffset = CS.DACOffset.get();
        scan.data.CSGain = CS.DACGain.get();
        scan.data.CSPreDelay=CS.PreDelay.get();
        scan.data.CSPostDelay=CS.PostDelay.get();
        scan.data.CSNSamps=CS.NSamps.get();
        scan.data.CSSampleCount=CS.SampleCount.get();
       d=smrun(real_scan,smnext('Ramsey_CDS_off_R'));       
       sm_getgradient
       if any(isnan(d{1})); break; end;
       d3(j,:)=nanmean(d{1});
   else
       d3(j,:) = NaN(1,npls);
      fprintf('no lock grad \n'); 
   end
    
end


%% make new CDS, RUN WITH CARE!!!

clear CS;
b=instrfind;
delete(b(12));
CS=CDS('COM3');
