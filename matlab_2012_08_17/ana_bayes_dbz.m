function out = ana_bayes_dbz(file, config)
%function out = ana_rescale_dbz_multi(file, config)
% analyzes many files taken in a row for dbz phase estimation
% uses a fourier transform to estimate dbz
% assumes all the files have the same mean dbz, fourier times, etc
dbzall=[]; dbz2all=[]; figure(66); clf; 
if ~exist('file','var') || isempty(file)
   file = get_files('sm*.mat'); 
end
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
%config = def(config,'fourier_inds',1:200); found later
config = def(config,'threshhold',0.4);
config = def(config,'nfreqs',250); 
config = def(config,'nmeas',2); 
% amass all of the data before processing
results = [];
for j = 1:length(file)
    s=load(file{j});
    data=s.data;
    scan=s.scan;
    scantime=getscantime(scan,data);
    if isfield(config,'t1')
       t1 = config.t1; 
    else
        [t1t t1] = att1(config.side,scantime,'after'); % FIX ME
    end
    data_all=anaHistScale(scan,data,t1); % FIX ME
    data=squeeze(data_all{1});
    xv = plsinfo('xval',scan.data.pulsegroups.name,[],scantime);
    if j ==1
        if isfield(scan.data,'pre_dbz') && isfield(scan.data,'post_dbz')
            config = def(config,'m_dbz',.5*1e-3*abs(scan.data.pre_dbz+scan.data.post_dbz)); %in GHz
            if abs(abs(scan.data.pre_dbz)-abs(scan.data.post_dbz)) > 15
               warning('gradient likely ran away!!!!') 
            end
        else
            if ~isfield(config,'m_dbz') %needed to stop asking for dbz
              config = def(config,'m_dbz',input('please enter m_dbz (in GHz)'));
            end
        end        
         config = def(config,'fourier_inds',1:(1+find(diff(diff(xv)))));
         %This only works for one or two repetitions. 
         if isopt(config,'rand') 
             config.fourier_inds=1:(length(xv)-50);          
         end
         if length(xv)>2*length(config.fourier_inds) && any(xv(config.fourier_inds)==xv(config.fourier_inds+length(config.fourier_inds)))             
             inds=find(xv(config.fourier_inds)==xv(config.fourier_inds+length(config.fourier_inds)));
             fourier_inds = [config.fourier_inds config.fourier_inds(inds)+length(config.fourier_inds)];                       
             config.fourier_inds=config.fourier_inds(inds)+length(config.fourier_inds).*(config.nmeas-1); 
             config.fourier_inds2=setdiff(fourier_inds,config.fourier_inds);
             config = def(config,'ftime2',xv(config.fourier_inds2))             
         else
             fourier_inds=config.fourier_inds; 
         end

         config = def(config,'ftime',xv(config.fourier_inds));
         config = def(config,'evo_inds',setdiff(1:length(xv),fourier_inds));
    end
       
    
    results(j).file = file{j};
    results(j).fourier_data = data(:,config.fourier_inds);    
    results(j).evo_data = data(:,config.evo_inds);
    results(j).t = xv(config.evo_inds);
    if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2)
        results(j).fourier_data2 = data(:,config.fourier_inds2);
    end 
end
config = def(config,'N_avg',500/2);N = size(s.data{1},1)/round(size(s.data{1},1)/config.N_avg);
D=size(results(1).evo_data,1)*size(results(1).evo_data,2)/N;
averager=kron(eye(D),ones(1,N))/N; %will average N sequential data points, matrix has D rows, ND columns

results2 = results;
tnow=tic;
avg_err = [];

if isopt(config,'calib')
    fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^2);
    beta0 = [.4, .35, 11, -pi/2, 600 0];
    evodata=(results(end).evo_data>config.threshhold);
    pars = fitwrap('plfit plinit',results(end).t,nanmean(evodata),beta0,fitfn, [1 1 1 1 1 0]);%fit oscillations
    dbz_avg=1./pars(3)
    calib(1)=1-2*pars(1);
    calib(2)=2*pars(2);
%     alpha=linspace(calib(1)-.3,calib(1)-.3,10);
%     beta=linspace(calib(2)-.8, calib(2)-.4,10);
%     for i=1:10
%         for j=1:10
%                 [dbz stddbz] = get_dbzbayes_act(results2(8).fourier_data,config.ftime,2*pi*0.012,2*pi*config.m_dbz,config.nfreqs,[alpha(i) beta(j)]);
%                 stdmat(i,j)=stddbz;
%                 dbzmat(i,j)=mean(dbz); 
%         end
%     end
%     stdmat
else
    calib=[0.25 0.67]; %calib(1)=a calib(2)=b
end

for j = 1:length(results)
    if isopt(config,'bayes')             
        T_samp = abs(diff(config.ftime(1:2))); F_samp = 1/T_samp;
        nyq = abs(F_samp/2); %nyquist frequency
        n_alias = floor(config.m_dbz/nyq); % number of times signal is aliased          
        [dbz stddbz] = get_dbzbayes_act(results2(j).fourier_data,config.ftime,2*pi./(2*T_samp),2*pi*config.m_dbz,config.nfreqs,calib);                
            if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2)            
                dbz2 = get_dbzbayes_act(results2(j).fourier_data2,config.ftime2,2*pi./(2*T_samp),2*pi*config.m_dbz,config.nfreqs,calib);                
                dbz2all=[dbz2all; dbz2];      
                 figure(41); hold on; 
                plot(dbz,dbz2,'.')
                dbzall=[dbzall; dbz];
            end
    elseif isopt(config,'bayslw')  
        T_samp = abs(diff(config.ftime(1:2))); F_samp = 1/T_samp;
        nyq = abs(F_samp/2); %nyquist frequency
        n_alias = floor(config.m_dbz/nyq); % number of times signal is aliased          
        dbz = get_dbzbayes(results2(j).fourier_data,config.ftime,F_samp,n_alias);
    else
        dbz = get_dbzthreshold(results2(j).fourier_data,config.ftime, config.m_dbz);
        dbzall=[dbzall; dbz]; 
        if isfield(config,'fourier_inds2') && ~isempty(config.fourier_inds2)
            dbz2 = get_dbzthreshold(results2(j).fourier_data2,config.ftime, config.m_dbz);    
            dbz2all=[dbz2all; dbz2];         
        end
        figure(41); hold on; 
        plot(dbz,dbz2,'.')
        
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
    %           aliased)
    
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

function [meandbz stddbz] = get_dbzbayes_act(fdata,ftime,B,dbz_avg,nfreqs,calib)
%a and b are correction terms. 

L = length(ftime); 
nalias=floor(dbz_avg/B);
ffreqs =(nalias+0.5).*B+B./2*linspace(-0.75,0.75,nfreqs); %choose frequencies to use. 

b=calib(2); %0.67; %correction terms. 
a=calib(1); %0.25; 
[wmat tmat]=meshgrid(ffreqs,ftime);
bayscoef=a+b*cos(wmat.*tmat); %update the bayesian. 

fdata1=2*(fdata>0.4)-1;
fdata2=repmat(fdata1,[1,1,length(ffreqs)]);
bayscoef2(1,:,:)=bayscoef;
bayscoef3=repmat(bayscoef2,[size(fdata2,1) 1 1]); %reshape to work w/ fdata3;
baysup=1-bayscoef3.*fdata2;

postdist=squeeze(prod(baysup,2));
postdist=bsxfun(@rdivide,postdist,sum(postdist,2));
ffreqsrep=repmat(ffreqs,[size(fdata,1),1]);
meandbz=sum(ffreqsrep.*postdist,2);

meandbz=meandbz./(2*pi); 
stddbz=nanstd(meandbz); 

figure(66); hold on; 


num=floor(size(postdist,1).*rand);
pd=postdist'; 
plot(ffreqs',pd(:,num)); 

end