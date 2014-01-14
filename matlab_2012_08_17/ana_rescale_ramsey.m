function out = ana_rescale_ramsey(file, config)
%function out = ana_rescale_ramsey(file, config)
% analyzes files with some period of frequency estimation, either dbz or
% ramsey, and then uses these frequencies to rescale ramsey fringes.
% 
% Uses Fourier transform or baysian estimator to estimate freqs, assumes 
% all the files have the same mean dbz, fourier times, etc
% 
% Assumes the pulsegroup used the option 'multi' and looks for two varpars
% to distinguish between the estimation and evolution segments.

% Possible opts:
% 'ramsey_ramsey': assumes ramsey data were taken for estimation. Default
% behavior assumes that dbz data were taken for estimation.
% 'bayes' uses baysian algorithm to compute frequencies
% 'fpga' uses fpga data to in place of alazar data
% 'nofit' does not attempt to fit the data, useful for finding the right
% detuning range for dbz estimation.
% 
% Potentially non-obvious config fields
% 'fpga_thresh': overrides the fpga threshold saved with the data set.
% 'fit_range': Range of times to fit. Useful to exclude outliers afte
% reconstruction.
% 'N_avg': number of points to average together after data has been
% rescaled.

dbzall=[]; dbz2all=[]; figure(66); clf; 
if ~exist('file','var') || isempty(file)
   file = get_files('sm*.mat'); 
end

%Flip the files to make sure the shortest evo time is first.
file=fliplr(file);

if isempty(file) || (isnumeric(file{1}) && file{1}==0)
    return
end

if  ~exist('config','var')
    config = struct();
end
figure(41); clf; 

out.file = file;
config = def(config,'opts','');
config = def(config,'decay',1.5);
config = def(config,'rng',[1 inf]);
config = def(config,'side','right');
config = def(config,'figind',470); figind = config.figind; fb = figind;
config = def(config,'threshhold',0.4);
config = def(config,'nfreqs',512); 
config = def(config,'nmeas',2); 
config = def(config,'threshold',0.4); 
config = def(config,'startind',2); %Default is 2 because need an extra pulse to trigger CDS
config = def(config,'fit_range',[100 2000]); % Only include times after this in the fit
config = def(config,'phse',pi); 
config = def(config,'burn',0); 

% amass all of the data before processing
results = [];
for j = 1:length(file)
    s=load(file{j});
    data=s.data; scan=s.scan;
    scantime=getscantime(scan,data);
    if isfield(config,'t1')
       t1 = config.t1; 
    else
        [t1t t1] = att1(config.side,scantime,'after'); %find t1 for histogramming. 
    end
    data_all=anaHistScale(scan,data,t1); % Histogram all data; 
    data=squeeze(data_all{1});
    xv = plsinfo('xval',scan.data.pulsegroups.name,[],scantime);
    if j ==1
        if isfield(scan.data,'pre_dbz') && isfield(scan.data,'post_dbz') %set the mean nuclear gradient based on first scan. 
            %config = def(config,'m_dbz',.5*1e-3*abs(scan.data.pre_dbz+scan.data.post_dbz)); %in GHz
            config = def(config,'m_dbz',1*1e-3*abs(scan.data.post_dbz)); %in GHz
            if abs(abs(scan.data.pre_dbz)-abs(scan.data.post_dbz)) > 15
               warning('gradient likely ran away!!!!') 
            end
        else
            if ~isfield(config,'m_dbz') %needed to stop asking for dbz
              config = def(config,'m_dbz',input('please enter m_dbz (in GHz)'));
            end
        end  
        
        if isopt(config,'ramsey_ramsey')
            config.m_dbz=0.005;
        end
        
        %only works with pulsegroups with two pulses.
        gd=plsinfo('gd',s.scan.data.pulsegroups(1).name,[],scantime);
        fourier_size=size(gd.varpar{1},1);
        evo_size=size(gd.varpar{2},1);
        
         startind=config.startind;
         %config = def(config,'fourier_inds',startind:(1+find(diff(diff(xv))))); 
         config = def(config,'fourier_inds',startind:fourier_size);
         config = def(config,'evo_inds',fourier_size+1+config.burn:fourier_size+evo_size); 
         %This only works for one or two repetitions. 
         if isopt(config,'rand') 
             config.fourier_inds=1:(length(xv)-50);          
         end
%          if length(xv)>2*length(config.fourier_inds) && any(xv(config.fourier_inds)==xv(config.fourier_inds+length(config.fourier_inds)))             
%              inds=find(xv(config.fourier_inds)==xv(config.fourier_inds+length(config.fourier_inds)));
%              fourier_inds = [config.fourier_inds config.fourier_inds(inds)+length(config.fourier_inds)];                       
%              config.fourier_inds=config.fourier_inds(inds)+length(config.fourier_inds).*(config.nmeas-1); 
%              config.fourier_inds2=setdiff(fourier_inds,config.fourier_inds);
%              config = def(config,'ftime2',xv(config.fourier_inds2))             
%          else
%          end
         
         fourier_inds=config.fourier_inds; 

         config = def(config,'ftime',xv(config.fourier_inds));
%          config = def(config,'evo_inds',setdiff(startind:length(xv)-(startind-1),fourier_inds));
    end
       
    
    results(j).file = file{j};
    results(j).fourier_data = data(:,config.fourier_inds);    
    results(j).evo_data = data(:,config.evo_inds);
    results(j).t = xv(config.evo_inds);
    if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2)
        results(j).fourier_data2 = data(:,config.fourier_inds2);
    end 
    if isfield(scan, 'data') && isfield(scan.data,'FPGA') && isfield(scan.data.FPGA,'data') && ~isempty(scan.data.FPGA)
        if any(fliplr(size(scan.data.FPGA.data))==size(results(j).fourier_data))
            results(j).fpga_data=scan.data.FPGA.data;
        else
            results(j).fpga_data=scan.data.FPGA.data';
        end
    end
    
    if isopt(config,'fpga') && any(~(size(results(j).fpga_data')==size(results(j).fourier_data)))
        error('Fourier data and FPGA data are not the same size');
    end
    if isopt(config,'fpga') %Use FPGA data in place of Fourier data.
        if isfield(scan.data,'FPGAThreshold') && ~isempty(scan.data.FPGAThreshold)
            config=def(config,'fpga_thresh',scan.data.FPGAThreshold);
        elseif isfield(scan.data.FPGA,'Threshold') && ~isempty(scan.data.FPGA.Threshold)
            config=def(config,'fpga_thresh',scan.data.FPGA.Threshold);
        end
        results(j).fourier_data = results(j).fpga_data'<config.fpga_thresh;
    end
     if isopt(config,'fpgaold') %Histogram the fourier_data and the evo_data separately.
        data_fourier=s.data; 
        data_fourier{1}=-1*data_fourier{1}(:,config.fourier_inds); 
        data_fourier=anaHistScale(scan,data_fourier,t1); % FIX ME
        results(j).fourier_data=squeeze(data_fourier{1});
        
        data_evo=s.data; 
        data_evo{1}=data_evo{1}(:,config.evo_inds); 
        data_evo=anaHistScale(scan,data_evo,t1);
        results(j).evo_data=squeeze(data_evo{1});

    end
    if isopt(config,'noscaleold') %Manually threshold the fourier_data at the specified threshold, and auto scale evo_data
        results(j).fourier_data=squeeze(s.data{1}(:,config.fourier_inds));        
        results(j).fourier_data = results(j).fourier_data<config.fpga_thresh;
 
        data_evo=s.data; 
        data_evo{1}=squeeze(data_evo{1}(:,config.evo_inds)); 
        data_evo=anaHistScale(scan,data_evo,t1);
        results(j).evo_data=squeeze(data_evo{1});
    end
    if isfield(scan.data,'FPGA') && isfield(scan.data.FPGA,'freqs') && ~isempty(scan.data.FPGA.freqs)        
        results(j).FPGAfreqs=scan.data.FPGA.freqs;
    end
    
end
out.t1=t1;
results2 = results;
tnow=tic;
avg_err = [];

if isopt(config,'calib')
    fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^2);
    beta0 = [.25, .25, 11, -pi/2, 60 0];
    evodata=(results(1).evo_data>config.threshold);
    pars = fitwrap('plfit plinit',results(1).t,nanmean(evodata),beta0,fitfn, [1 1 1 1 1 0]);%fit oscillations
    dbz_avg=1./pars(3);
    calib(1)=1-2*pars(1);
    calib(2)=2*pars(2);
    fprintf('alpha is %f, beta is %f \n',calib(1),calib(2)); 

else
    calib=[0.25 0.67]; %calib(1)=a calib(2)=b
end

config = def(config,'T_samp',abs(diff(config.ftime(1:2)))); 
T_samp=config.T_samp; F_samp = 1/T_samp;
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(config.m_dbz/nyq); % number of times signal is aliased

%Estimate the frequencies from the estimation period.
for j = 1:length(results)
    if isfield(results,'FPGAfreqs')
        ff=@(x) (n_alias)*nyq+(-1)^(n_alias)*x/512*nyq+mod(n_alias,2)*nyq;
        results(j).FPGAfreqs=ff(results(j).FPGAfreqs);            
    end
    if isopt(config,'bayes')        
        config = def(config,'bayes_freq',[nyq*n_alias,nyq*(n_alias+1)]);         
        bayes_freq=config.bayes_freq; 
        if bayes_freq(1)>0.3 
            warning('This is a large frequency. Please input frequency in GHz')
        end
        [dbz stddbz,stddbzall,badinds] = get_dbzbayes_act(results2(j).fourier_data,config.ftime,bayes_freq,config.nfreqs,calib,config.threshold,config.phse);                
        dbzall=[dbzall; dbz]; %dbzall and dbz2all collect all the mean values measured. 
        results(j).badinds=badinds;
        %returns mean and std dev. of dbz
        %args: data, times of sampling, range of distribution, center of distribution, # freq in distribution ,alpha/beta vals     
            if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2) %if there are two measurement sets in each pulse, this analyzes those            
                dbz2 = get_dbzbayes_act(results2(j).fourier_data2,config.ftime2,pi./(4*T_samp),2*pi*config.m_dbz,config.nfreqs,calib);                
                figure(41); hold on; 
                plot(dbz,dbz2,'.')                 
                dbz2all=[dbz2all; dbz2]; 
            end
    elseif isopt(config,'bayslw')  
        dbz = get_dbzbayes(results2(j).fourier_data,config.ftime,F_samp,n_alias);
    elseif isopt(config, 'raw') %plot the raw data only, no correction
        dbz=config.m_dbz*linspace(1,1,size(results2(j).fourier_data,1))';      
    elseif isopt(config,'fpgafreqs') 
        dbz=results(j).FPGAfreqs;   
    else
        dbz = get_dbzthreshold(results2(j).fourier_data,config.ftime, config.m_dbz);
        dbzall=[dbzall; dbz]; 
        if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2)
            dbz2 = get_dbzthreshold(results2(j).fourier_data2,config.ftime, config.m_dbz);    
            dbz2all=[dbz2all; dbz2];         
            figure(41); hold on; 
            plot(dbz,dbz2,'.')
        end
    end

   
    results(j).dbzs=dbz;
    

end

config = def(config,'N_avg',size(results(1).evo_data,1)); 
%config = def(config,'N_avg',256); 

%Sweep through VCO freqs, and reconstruct for all data sets separately and together. Store all of the
%data for later analysis.

clear tau;
clear period;
clear const;
clear amp;

%default assumes dbz_ramsey, i.e., dbz estimation and ramsey evo.
if isopt(config, 'ramsey_ramsey')
    vcofreqs=0;  
else
    vstart=.0625;
    vend=.0675;
    vpoints=50;
    vcofreqs=linspace(vstart,vend,vpoints);
    vcofreqs(end+1)=vcofreqs(35);
end

for k=1:length(vcofreqs);
    
    %Find the mean detuning across all data sets
    for j = 1:length(results)
        
        if isopt(config, 'ramsey_ramsey')
            results(j).dbzs_ramsey=results(j).dbzs;
        else
            results(j).dbzs_ramsey=abs(results(j).dbzs-vcofreqs(k));
        end
   
    end
    
    evobig=[];
    tbig=[];
    
    %Rescale the data sets individually and together
    for j=1:length(results)
        
        %when doing post selection, will start by not changing mdbz. let's consider
        %this though...
        mdbz=mean(mean([results.dbzs_ramsey]));
        
    %Rescale and average each dataset individually
    evodata=(results(j).evo_data>config.threshhold); %threshold from histogramqubitdetection
    tstemp=results(j).dbzs_ramsey*results(j).t/mdbz;%this is a matrix of corrected times, with each entry corresponding to an entry of the evolution data, it's a column vector of magnetic field corrections for each of the 1000 runs multiplied by a row vector of the 50 nominal times
    if isopt(config,'postsel')
        evodata(results(j).badinds,:)=[];
        tstemp(results(j).badinds,:)=[];
    end

    tst=tstemp(:);%arrange all data points into a vector, the vector will be roughly time ordered
    evost= evodata(:);%similarly arrange the readout data, the entries of this vector will still correspond to the times in tst
    [tt ind]=sort(tst); %find correct time order
    avg_err = [avg_err,std(diff(tt))];
    evofinal=evost(ind); %sort measurement results into time order
    
    D=size(evofinal,1);
    D_new=floor(D/config.N_avg);
    
    % simple for-loop averager
    for m=1:D_new
        results(j).newdata(k,m) = nanmean(evofinal((m-1)*config.N_avg+1:m*config.N_avg));%average results into groups of N and then put into data structure
        results(j).newtime(k,m)=nanmean(tt((m-1)*config.N_avg+1:m*config.N_avg));
    end
    
    evobig=cat(1,evost,evobig);
    tbig=cat(1,tst,tbig);
    [tbigfinal ind]=sort(tbig); 
    evobigfinal=evobig(ind);
    
    %Rescale and average all datasets together.
    D=size(evobigfinal,1);
    D_new=floor(D/config.N_avg);
    
    % simple for-loop averager
    for m=1:D_new
        out.newdata(k,m) = nanmean(evobigfinal((m-1)*config.N_avg+1:m*config.N_avg));%average results into groups of N and then put into data structure
        out.newtime(k,m)=nanmean(tbigfinal((m-1)*config.N_avg+1:m*config.N_avg));
    end
    
    end
    



out.results=results; 
out.config = config;
out.detunings=vcofreqs;
boundary=round(.05*size(out.newtime(1,:),2)); 
fit_range(1) = size(out.newdata(1,:),2)-size(out.newtime(1,out.newtime(1,:)>config.fit_range(1)),2);
fit_range(2) = size(out.newtime(1,out.newtime(1,:)<config.fit_range(2)),2);

if isopt(config,'multi') && ~isopt(config, 'nofit')
    clear tau;
    clear const;
    clear period;
    allpars=zeros(50,length(results)+6);  
    for k=1:length(vcofreqs)
        data = [];
        model = [];
        fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
        mask = [1 1 1 1 0 0];
        beta0 = [.36, .36, 150, 500, 0, 1.5];
        boundary = 0;
        for j = 1:length(results)
            if all((~isnan(results(j).newdata)))
                mask2 = 3:(size(out.results(j).newdata,2)-2);                
                data(end+1).y = out.results(j).newdata(k,mask2);
                data(end).x = out.results(j).newtime(k,mask2);
                model(end+1).fn = fitfn;
                model(end).pt = @(p) [p(1:6), p(7+boundary)];
                boundary = boundary+1;
                mask = [mask 1];
                beta0 = [beta0, .01];
            end
        end
        
        pars = mfitwrap(data,model,beta0,'plinit plfit lm',mask);
        mask(6) = 1;
        pars = mfitwrap(data,model,pars,'plinit plfit lm',mask);
        allpars(k,:)=pars; 

        tau(k)=abs(pars(4));
        const(k)=pars(6);
        period(k)=pars(3);

    end
    figure(3); clf; plot(tau);
    xlabel('Detuning index');
    ylabel('T2* (ns)');
    [~,ind]=max(tau);
    figure(figind); figind = figind+1; clf; hold on;
    %for j = 1:length(data)
    for j=1:3
        plot(data(j).x,out.results(j).newdata(ind,mask2),'x-'); hold on; 
        plot(data(j).x,fitfn(model(j).pt(allpars(k,:)),data(j).x),'k');
    end
    xlabel('Evolution Time (ns)');
    ylabel('P_{Triplet}');
    title(sprintf('T2* is %f', tau))

elseif ~isopt(config,'nofit')
    %setup the mfit
    data = [];
    model = [];
        
    fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(6))).*exp(-((x)/p(4)).^p(5));
    pars = [.3, .3, 115, 1500, 2,-pi/2];
    lb = [.01 .01 75 50 1.5  0];
    ub=  [.7 .7 115 5000 3 2*pi];
        

    mask2 = max([1*boundary fit_range(1)]):min([fit_range(2) (size(out.newdata,2)-1*boundary)]);
    data(end+1).y = out.newdata(k,mask2);
    data(end).x = out.newtime(k,mask2);
    options = optimset('Display','off','MaxIter',10000,'TolFun',1e-13,'Algorithm',{'levenberg-marquardt',.005});
    %pars=lsqcurvefit(fitfn,pars,data(1).x,data(1).y,lb,ub,options);
    pars=lsqcurvefit(fitfn,pars,data(1).x,data(1).y,[],[],options);
    out.pars=pars; 
    fit=fitfn(pars,data(1).x);
    tau(k)=abs(pars(4));
    const(k)=pars(5);
    period(k)=pars(3);
    amp(k)=pars(2);
    
    if isopt(config, 'ramsey_ramsey')
    figure(1); clf; 
    plot(data(1).x,data(1).y, '.-', data(1).x, fit, 'k')
    titlestring=sprintf('T2* = %.5f, Decay const = %.5f',tau(k),const(k));
    title(titlestring);
    else
    figure(1); clf; hold on;
    subplot(4,1,1);
    plot(data(1).x,data(1).y, '.-', data(1).x, fit, 'k')
    xlabel('VCO index');
    ylabel('Triplet probability');
    subplot(4,1,2); plot(tau);
    xlabel('VCO index');
    ylabel('T2* (ns)');
    subplot(4,1,3); plot(amp);
    xlabel('VCO index');
    ylabel('Fit amplitude');
    subplot(4,1,4); plot(const);
    xlabel('Detuning index');
    ylabel('Decay constant');
    end

end
end
out.fit.tau=tau;
out.fit.const=const;
out.fit.period=period;

%plot a histogram of the frequencies fit.
dbzset=[results.dbzs];
[ct,bins]=hist(dbzset(:),30);
figure(44); clf; subplot(1,2,1); 
plot(bins,ct,'.-'); hold on;
title('Histogram of Frequencies measured')
legend('Computer','FPGA')
if isfield(results,'FPGAfreqs')
    fpgafreqs=[results.FPGAfreqs]; 
    [ct,bins]=hist(fpgafreqs(:),30);
    plot(bins,ct,'r.-')
    subplot(1,2,2); 
    plot(dbzset(:),fpgafreqs(:),'.') 
end



end


function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
end 

function pars = fitdbzoscillations(x,y,fitopts)
%cosfn = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x))';
cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
%cosfn3 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)) * y(6))';
cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
mask = [1 1 1 1 0 1];
fifn.fn = @fioscill;
fifn.args = {2};
fp=fitwrap('plinit plfit woff',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
if ~isempty(strfind(fitopts,'badinit'))
    pars = fitwrap(fitopts,x,y,[.2 .2 .2 .06 0 .01],cosfn2,mask);
else
   pars=fitwrap(fitopts,x,y,fp, cosfn2, mask);
end
end

function pars = fitgauss(x,y,fitopts)
fitfn = @(p,x) p(1)*exp(-(x-p(2)).^2/p(3)^2);
[~,mi]=max(smooth(y(5:end-5)));
beta0 = [.1*max(y)/range(x), x(mi),.1*range(x)];
pars = fitwrap(fitopts,x,y,beta0,fitfn);
end

function ff = getfreqs(Fsamp, NFT,mdbz)
nyq = abs(Fsamp/2);
n_alias = floor(mdbz/nyq);
ff = (n_alias+.5)*nyq+(-1)^n_alias*(Fsamp/2)*linspace(-.5,.5,NFT/2+1);
end

function s=def(s,f,v)
  if(~isfield(s,f))
     s.(f) = v;
  end
end

function dbzs = get_dbzthreshold(fdata,times,mdbz)
%function dbzs = get_dbz(fdata,times,mdbz)
    % estimates the frequency of dbz based on fourier data. inputs are
    % fdata: the fourier data: N_reps x number of different evo times used 
    %           N_reps different frequencies will be estimated
    % times: the different evolution times used to estimate it
    % mdbz: the mean value of dbz (used to figure out how many times
    %           aliased)    
    %version uses thresholded data
    
L = length(times);
T_samp = abs(diff(times(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
dbzs=zeros(size(fdata,1),1);
fdata=(fdata>0.4); %threshold comes from histogram analysis in histogramqubitdetection.m
for j = 1:size(fdata,1)
   ft = fft(fdata(j,:)-mean(fdata(j,:)),NFFT)/L;
   ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbzs(j) = ff(mi+1); %+1 because looking at 2:end
end

end

function meandbz = get_dbzbayes(fdata,times,B,nalias)
%function dbzs = get_dbzbayes(fdata,times,mdbz)
    % estimates the frequency of dbz based on fourier data. inputs are
    % fdata: the fourier data: N_reps x number of different evo times used 
    %           N_reps different frequencies will be estimated
    % times: the different evolution times used to estimate it
    % B: bandwidth of estimation (used to figure out how many times   
    %version uses thresholded data
    %version implements Bayes inversion
    
    %assume that data is measuring each time once


L = length(times);
K=L*(L+1)/2; %K from Steve's notes, size of Fourier Series vector
Kp=2*K+1;
% T_samp = abs(diff(times(1:2))); F_samp = 1/T_samp;
% %NFFT = 2^nextpow2(8*L); %padded with zeros
% nyq = 0.05; %nyquist frequency


freqs=-K:K;

%parameters for Bayes model
b=0.67;
%b=0.73; %values from overall data
a=0.25;

%setup to find mean and variance of posterior distribution
meanvector=(B*((-1).^freqs-1)./(freqs*pi).^2); % averaging for nonzero frequencies from equation 18, multiplied by 1/B which is part of the normalisation of pq
meanvector(K+1)=0; %division by zero turns up NaN that needs fixing


meanvector=(B*((-1).^freqs-1)./(freqs*pi).^2); % averaging for nonzero frequencies from equation 18, multiplied by 1/B which is part of the normalisation of pq
meanvector(K+1)=0; %division by zero turns up NaN that needs fixing


meandbz=zeros(size(fdata,1),1);
fdata=2*(fdata<0.4)-1; %threshold comes from histogram analysis in histogramqubitdetection.m, output is now ones and -1's, 1 corresponds to a Singlet outcome which is consistent with convention in our notes


R=size(fdata,1);

%update posterior distribution using Baye's rule
parfor j=1:R

    pq=zeros(Kp,1); %new posterior distribution
    pq(K+1)=1; % initially uniform
    for k=1:L
        pqold=pq;
        %pq=pqold+fdata(l,k)*(a*pqold+0.5*b*spdiags(ones(Kp,2),[-k, k],Kp,Kp)*pqold); 
        pq=pqold+fdata(j,k)*(a*pqold+0.5*b*(pqold([(k+1):Kp, 1:k])+pqold([(Kp-k+1):Kp, 1:(Kp-k) ]))); %implements update rule in notes equation 17, avoids matrix multiplication just adds vectorsfigure
        %there may be a better way to do this since i am assuming that
        %pq(1:k) and pq(Kp-k+1:Kp) are zero vectors, this is satisfied in
        %practice since our pq vector is bigger than it needs to be for
        %some reason. this method runs an order of magnitude faster than
        %the previous one, it may be better to build vectors with explicit
        %zeros but that may mess with the timing
    end
    meandbz(j)=((nalias+0.5)*B-meanvector*pq/pq(K+1))/2/pi; %division by DC component pq(K+1) handles normalisation of pq, factor of 1/B already accounted for, 2pi is there because our note is angular frequency everywhere
%minus sign here result of current fudge to reflect the answer around
%within its bandwidth, this is to do with some error in the notes
end
end

function [meandbz,stddbz,stddbzall,badinds] = get_dbzbayes_act(fdata,ftime,freq,nfreqs,calib,thresh,phse)
%The Bayesian updater has form P(w|r)=1-r(a +b cos (w t)) %
%Since a+b cos (wt) does not depend on results, create that set of n_w x
%n_t data points first (bayscoef). 
%Then threshold data so that T->1, S-> -1 (fdata). 
%wnat prod(1-fdata*bayscoef), but to get vectorizing right, must do a bunch
%of repmat/permute typestuff. 
%finally, finds freq. w/ max probability. 
%Plots histogram of standard deviations in MHz. 

%ffreqs =2*pi*nyq*linspace(n_alias,n_alias+1,nfreqs); %choose frequencies to use. 
ffreqs =2*pi*linspace(freq(1),freq(2),nfreqs); %choose frequencies to use. 
%ffreqs=2*pi*nyq*linspace(-.125,.125,nfreqs)+2*pi*.063;

b=calib(2); %0.67; %correction terms. 
a=calib(1); %0.25; 
[wmat,tmat]=meshgrid(ffreqs,ftime);
bayscoef=(a+b*cos(wmat.*tmat+phse));%.*exp(-tmat/10); %update the bayesian.  %this is n_t (rows) x n_w (cols)

fdata1=(2*(fdata>thresh)-1); %size nreps x nt
fdata2=repmat(fdata1,[1,1,length(ffreqs)]); %repeat the measurement data over the set of frequencies on 3rd dim. 
bayscoef2(1,:,:)=bayscoef; %dim 1 x n_t x n_w 
bayscoef3=repmat(bayscoef2,[size(fdata2,1) 1 1]); %reshape to work w/ fdata2, have time/freq array repeat along first dimensions;
baysup=1-bayscoef3.*fdata2; %bayes updater, multiplies variables (nreps x n_tms x dummy)*(dummy x n_tms x n_freqs). 

postdist=squeeze(prod(baysup,2)); %take the product along meas time dimension, end up with (nreps x nfreqs)
postdist=bsxfun(@rdivide,postdist,sum(postdist,2)); %normalize each distribution by taking sum along freq. 
postdistcum=cumprod(baysup,2); 
postdistcum=bsxfun(@rdivide,postdistcum,sum(postdistcum,3));
[~,max_dbz_time]=max(postdistcum,[],3); 
%mval=mean(max_dbz_time(:,end-10:end),2);
mval=max_dbz_time(:,end-10); 


ffreqsrep=repmat(ffreqs,[size(fdata,1),1]); %dummy for nreps x nfreqs w/ freq vals)
%meandbz=sum(ffreqsrep.*postdist,2);
[~,meandbz_ind]=max(postdist,[],2); 
diffind=abs(mval-meandbz_ind);
badinds=find(diffind>2); 
meandbz=ffreqs(meandbz_ind)'; 
stddbzall=(sum(ffreqsrep.^2.*postdist,2)-sum(ffreqsrep.*postdist,2).^2).^(0.5); 

meandbz=meandbz./(2*pi); 

stddbz=nanstd(meandbz); 

figure(66); hold on; 
[ct,bins]=hist(stddbzall,50);


plot(bins/(2*pi)*1000,ct); 
xlabel('Frequency (MHz)'); 
ylabel('Count')
%num=floor(size(postdist,1).*rand)+1;
pd=postdist'; 
%plot(ffreqs',pd(:,num)); 

end



