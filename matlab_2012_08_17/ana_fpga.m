%% Load data
fpgadata=importdata('Z:\\John\FPGA data\\fpga_array_0.dat');
% fpgadata=fpgadata'*-1;
% fpgafreqs=importdata('Z:\\John\FPGA data\\data3.dat');

% d=ana_avg();
% alazardata=d.data{1}(:,2:200);
% close all;

% figure(3); plot(alazardata(5,:)); hold on; plot(fpgadata(5,:),'r');

% for j=1:size(fpgadata,1)
%     C=corrcoef(fpgadata(j,:),alazardata(j,:));
%     cor(j)=C(1,2);
% end
% 
% % figure(1);
% % plot(cor);
% 
% % figure(2);
% % plot(fpgafreqs);
% 
% %Threshold
% alazardatat=-(2*(alazardata>.4)-1);
fpgadata=fpgadata';
fpgadatat=1*(2*(fpgadata>-3500)-1);

ftime = 12:12:199*12;

% dbza = get_dbzbayes_jmn_v2(alazardata,ftime,2*pi*0.05,1);
dbzf=get_dbzbayes_jmn_nothreshold(fpgadatat,ftime,2*pi*0.05,1);
dbzf2=get_dbzbayes_jmn_nothreshold_v2(fpgadatat,ftime,2*pi*0.05,1,.06);
%temp2=2* 4.1667e+07-dbzf* 4.1667e+07/256;
figure(5); clf;
plot(dbzf); 
%yaxis([.04 .09])
figure(7); clf;
plot(dbzf2); 
yaxis([.04 .09])


% figure(3); clf;
% plot(dbzba,'b'); hold on; plot(fpgafreqs,'r');

%% Run through data and look at histograms and raw data

for j=1:size(alazardata,1)
    figure(1);
    clf;
    hold on;
    
    %hist(alazardata(j,:),50);
%     subplot(2,1,1);
%     plot(alazardata(j,:),'.');
%     yaxis([-2 2]);
%     subplot(2,1,2);
%     plot(fpgadata(j,:),'.');
    ma(j)=mean(alazardata(j,:));
    mf(j)=mean(fpgadata(j,:));
    yaxis([2000 5000]);
%     pause(.5);
    
end

%% Repackage the FPGA data to look like the real data for analysis with ana_rescale etc.

temp = get_files('sm*.mat');
files=fliplr(temp);

for j=1:length(files);
  %fname=sprintf('Z:\\John\\FPGA data\\fpga_array_%d.dat',j-1);
  fname=sprintf('fpga_12_20_array_%d.dat',j-1);
  fpgadata=load(fname);
  fpgadata=fpgadata';
  
  d=ana_avg(files{j},'noscale');
  fname
  files{j}

%Fit alazar vs fpga to rescale  
%   data=squeeze(d.data{1});
%   alazardata=data(:,2:200);
%   p=polyfit(fpgadata(:),alazardata(:),1);
%   fpgadata=fpgadata*p(1)+p(2); 
%   
%   figure(51); clf;
%   hist(d.data{1}(:),50);

  d.data{1}(1:1024,1,2:200)=fpgadata(1:1024,1:199);
  
%   figure(52); clf;
%   hist(fpgadata(:),50);
  
  savename=sprintf('sm_raw_FPGA__12_20_%d.mat',j);
  d.filename=savename;
  save(savename,'-struct','d');
 
  
end

%%
    figure(1); clf;
    hold on;
for j=1:19;

    plot(aout.results(j).t,mean(aout.results(j).evo_data));
end

figure(2); clf
hold on;
for j=1:19;

    plot(fout.results(j).t,mean(fout.results(j).evo_data));
end

%%
fname=sprintf('Z:\\John\\FPGA data\\fpga_freqs_%d.dat',1-1);
        freqs=importdata(fname);
        actfreqs=2* (1/12/2)-freqs* (1/12/2)/256;
        figure(9); plot(actfreqs)
        
%%
d=ana_avg();
data=squeeze(nanmean(d.data{1}));
figure(25); clf;
plot(data(181:end));
title('G=-20, O=23000');

%%
d=ana_avg();
fourier_data=d.data{1}(:,2:200);
evo_data=d.data{1}(:,201:600);
evo_data=fliplr(evo_data);

%% Fit reconstructed ramsey fringes
clear tau;
clear const;
for k=25:75;
%setup the mfit
data = [];
model = [];
%Important to get p(7), which is the phase correct.
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
mask = [1 1 1 1 0 0 0];
beta0 = [.45, .3, 150, 500, 0, 1.1,pi/2];
offset = 0;
mask2 = 10:(size(out.results(1).newdata,2)-75);
data(end+1).y = out.results(1).newdata(k,mask2);
data(end).x = out.results(1).newtime(k,mask2);
% Fit the parameters gradually.
pars=fitwrap('',data(1).x,data(1).y,pars,fitfn,mask);
mask(7) = 1;
pars=fitwrap('',data(1).x,data(1).y,pars,fitfn,mask);
mask=0*mask;
mask(3)=0;
mask(6)=0;
mask(4)=1;
pars=fitwrap('plfit robust',data(1).x,data(1).y,pars,fitfn,mask);
tau(k)=abs(pars(4));
const(k)=pars(6);
figure(2); clf; plot(tau);
xlabel('Detuning index');
ylabel('T2* (ns)');
end

%% Fit reconstructed ramsey fringes using lsqcurvfit
%x = lsqcurvefit(fun,x0,xdata,ydata)
% x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)

clear tau;
clear const;
clear period;
for k=25:75;
%setup the mfit
data = [];
model = [];
%Important to get p(7), which is the phase correct.
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(6))).*exp(-((x)/p(4)).^p(5));

beta0 = [.45, .3, 150, 500, 0, 1.1,pi/2];
lb = [0 0 100 50 0  1.0 0];
ub=  [.5 .5 400 2000 1.6 2*pi];
mask2 = 10:(size(out.results(1).newdata,2)-50);
data(end+1).y = out.results(1).newdata(k,mask2);
data(end).x = out.results(1).newtime(k,mask2);
options = optimset('Display','off','MaxIter',1000);
pars=lsqcurvefit(fitfn,beta0,data(1).x,data(1).y,lb,ub,options);
fit=fitfn(pars,data(1).x);
figure(1); clf;
plot(data(1).x,data(1).y, '.-', data(1).x, fit, 'k')
tau(k)=abs(pars(4));
const(k)=pars(6);
period(k)=pars(3);
figure(3); clf; plot(tau);
xlabel('Detuning index');
ylabel('T2* (ns)');
end


%% Plot the reconstructed ramsey fringes
set=70;
figure(1); clf; plot(out.results(1).newtime(set,mask2),out.results(1).newdata(set,mask2));

%%
figure(2); clf;
index=25;
plot(out.newtime(index,:),out.newdata(index,:)) 

%%
for j=100:length(out.detunings)
    figure(1); clf; plot(out.newtime(j,:),out.newdata(j,:)) 
    titlestring=sprintf('%.5f',out.detunings(j) );
    title(titlestring);
    pause(.3);
end

%% Plot histograms from Alazar and FPGA

file = get_files('sm*.mat');
d=load(file{1});
figure(1); clf; hist(d.data{1}(:),50);
%figure(2); clf; hist(d.scan.data.FPGA.data(:),50);
figure(3); clf; imagesc(squeeze(d.data{1}));
figure(4); clf; imagesc(squeeze(d.data{1}));

%%
offset=round(.05*size(out.newtime(1,:),2));
figure(1); clf;
clear tau;
clear period;
clear const;
clear amp;
for k=1:length(out.detunings)
    %setup the mfit
    data = [];
    model = [];
        
    fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(6))).*exp(-((x)/p(4)).^p(5));
    pars = [.45, .3, 200, 500, 0, 1.6,pi/2];
    lb = [.1 .1 75 50 1.1  0];
    ub=  [.5 .5 500 2000 1.6 2*pi];
        
    mask2 = offset:(size(out.newdata,2)-offset);
    data(end+1).y = out.newdata(k,mask2);
    data(end).x = out.newtime(k,mask2);
    options = optimset('Display','off','MaxIter',1000);
    pars=lsqcurvefit(fitfn,pars,data(1).x,data(1).y,lb,ub,options);
    fit=fitfn(pars,data(1).x);
    tau(k)=abs(pars(4));
    const(k)=pars(5);
    period(k)=pars(3);
    amp(k)=pars(2);
    figure(1); clf; hold on;
    subplot(4,1,1);
    plot(data(1).x,data(1).y, '.-', data(1).x, fit, 'k')
    subplot(4,1,2); plot(tau);
    xlabel('Detuning index');
    ylabel('T2* (ns)');
    subplot(4,1,3); plot(amp);
    xlabel('Detuning index');
    ylabel('Amplitude');
    subplot(4,1,4); plot(const);
    xlabel('Detuning index');
    ylabel('Decay exponent');

end

%%
%% plot data sets with FPGA parameters
d = ana_avg;
d=flipud(d)
gain=[];
close all;
figure(668); clf; hold on; 
figure(666); clf; hold on; c = ['k' 'r' 'g' 'b' 'c' 'm' 'y'];
results = struct();
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^p(6));
clear T2err; clear T2rng; clear results;  
c='krgbcmy';
ydata=[]; xdata=[];
for j =1:length(d)
    gd=plsinfo('gd',d(j).scan.data.pulsegroups(1).name,[],d(j).scantime);
    evo_size=0; 
    evo_st=find(gd.pulses==72); 
    varlengths=cellfun(@length,gd.varpar); 
    fourier_size=sum(varlengths(1:evo_st-1)); 
    evo_size=sum(varlengths(evo_st:end)); 
    evoinds=fourier_size+[6:evo_size]; 
    %param(j)=d(j).scan.data.FPGA.Offset;
    param(j)=d(j).scan.data.FPGA.Gain;
    cnum=mod(j,7)+1;
    %names{j}=sprintf('Offset %d file %s',param(j), d(j).filename(4:end-4));
    names{j}=sprintf('Gain %d file %s',param(j), d(j).filename(4:end-4));
    results(j).data = squeeze(nanmean(d(j).data{1}(:,:,evoinds)))';
    results(j).xv = d(j).xv{1}(evoinds);
    figure(666);
    if j<8
    plot(results(j).xv,results(j).data,'.-','Color',c(cnum),'DisplayName',names{j})        
    else
    plot(results(j).xv,results(j).data,'*-','Color',c(cnum),'DisplayName',names{j})        
    end
    

end

legend show;


%% Fit the data sets
data = [];
model = [];
figure(667); clf; hold on;
period=[]; tau=[]; betas=[];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
mask = [1 1 1 1 0 0 0];
boundary = 0;
for j = 1:length(results)
    cnum=mod(j,7)+1;
    if all((~isnan(results(j).data)))
        %maskfit=6:length(results(j).xv);
        maskfit=1:length(results(j).xv); 
        data(j).y = results(j).data(:,maskfit); 
        data(j).x = results(j).xv(:,maskfit);
        fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(6))).*exp(-((x)/p(4)).^p(5));
        pars = [.3, .3, 2000, 1500, 2,-pi/2];
        lb = [.01 .01 75 50 1.5  0];
        ub=  [.7 .7 500 5000 3 2*pi];
        options = optimset('Display','off','MaxIter',10000,'TolFun',1e-13,'Algorithm',{'levenberg-marquardt',.005});
        %pars=lsqcurvefit(fitfn,pars,data(1).x,data(1).y,lb,ub,options);
        pars=lsqcurvefit(fitfn,pars,data(j).x,data(j).y,[],[],options);
   
        fit(j,:)=fitfn(pars,data(j).x);
        tau(j)=abs(pars(4));
        period(j)=pars(3);
        tau(j)=pars(4);
        
        figure(667); hold on;
        plot(results(j).xv,results(j).data,'.-','Color',c(cnum),'DisplayName',names{j})    
        plot(results(j).xv,fit(j,:),'-','Color',c(cnum)) 
    end
    
end

legend show;


figure(667); clf; hold on;
    for j=1:2; 
        %plot(data(j).x,out.results(j).data(ind,mask2),'x-'); hold on; 
        plot(data(j).x,fitfn(model(j).pt(betas(j,:)),data(j).x),'r');
    end









