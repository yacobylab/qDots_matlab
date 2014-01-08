%%
%%

files = get_files('sm*.mat');
%%

for k= 1:length(files)
d = ana_avg(files{k});
data = squeeze(d.data{1});
mdbz = .5*(d.scan.data.post_dbz+d.scan.data.pre_dbz);
ttimes = d.xv{1}(1:200);
L = length(ttimes);
T_samp = abs(diff(ttimes(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);

fdata = data(:,1:200);
fdata=(fdata>0.4); %threshold comes from histogram analysis in histogramqubitdetection.m
dbz1=zeros(size(fdata,1),1);
dbz2=zeros(size(fdata,1),1);
for j = 1:size(fdata,1)
   ft = fft(fdata(j,:)-mean(fdata(j,:)),NFFT)/L;
   ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbz1(j) = ff(mi+1); %+1 because looking at 2:end
end
fdata = data(:,end-199:end);
fdata=(fdata>0.4); %threshold comes from histogram analysis in histogramqubitdetection.m
for j = 1:size(fdata,1)
   ft = fft(fdata(j,:)-mean(fdata(j,:)),NFFT)/L;
   ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbz2(j) = ff(mi+1); %+1 because looking at 2:end
end
results(k).dbz1 = dbz1;
results(k).dbz2 = dbz2;
results(k).wait = length(d.xv{1})-200;
end

%%
for j =1:length(results)
   results(j).ddbz = results(j).dbz1-results(j).dbz2;
   mask = abs(results(j).ddbz)< 1e-3;
   results(j).m_ddbz = mean(results(j).ddbz);
   results(j).sig_ddbz = std(results(j).ddbz(mask));
   [results(j).h, results(j).b]=hist(results(j).ddbz(mask),100);
   
end

%% pack it up for sydney
results  = [];
for j = 1:length(files);
   d = ana_avg(files{j});
   results(j).data = squeeze(d.data{1});
   results(j).times = d.xv{1};
   results(j).fname = d.filename;
   results(j).plslength = 4e-6;
end
%%
tic
d = ana_avg('sm_dBz_phase_raw_R_8217');
d2=anaRawScale(anaRawUnpack(d.scan,d.data),1e-2);
toc % 11s
d2{1} = squeeze(d2{1});
d2{1} = reshape(d2{1},length(d.xv{1}),d.scan.data.conf.nloop,d.scan.data.conf.nrep);
% should be npls x nrep x nloop

%%
mdbz = 1e-3*abs(d.scan.data.post_dbz+d.scan.data.pre_dbz)/2;
ttimes = d.xv{1}(1:200);
L = length(ttimes);
T_samp = abs(diff(ttimes(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
dbzs = zeros(size(d2{1},2),size(d2{1},3));

for j = 1:size(d2{1},3)
    fdata = d2{1}(:,:,j)-repmat(mean(d2{1}(:,:,j)),size(d2{1},2),1);
    ft = fft(fdata,NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1,:));
    [~, mi] = max(ft(2:end,:));
     dbz1(j,:) = ff(mi+1); %+1 because looking at 2:end
    
end


%%
tic
%d = ana_avg; 
d = ana_avg('sm_dBz_phase_raw_R_8217');close all
d2=anaRawScale(anaRawUnpack(d.scan,d.data),1e-2);
toc % 11s
d2{1} = squeeze(d2{1});

%%

mdbz = 1e-3*abs(d.scan.data.post_dbz+d.scan.data.pre_dbz)/2;
ttimes = d.xv{1}(1:200);
L = length(ttimes);
T_samp = abs(diff(ttimes(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
dbzs = zeros(1,size(d2{1},2));

for j = 1:size(d2{1},2)
   fdata = d2{1}(:,j) > 0.4;
   ft = fft(fdata-mean(fdata),NFFT)/L;
   ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbzs(j) = ff(mi+1); %+1 because looking at 2:end

end
dbz2 = reshape(dbzs,d.scan.data.conf.nloop,1024);


for j = 1:size(dbz2,2)
   mask = abs(dbz2(:,j)-mean(dbz2(:,j)))< .002;
   s_dbz(j) = std(dbz2(mask,j));
end

%%
mdbz = 1e-3*abs(d.scan.data.post_dbz+d.scan.data.pre_dbz)/2;
ttimes = d.xv{1}(1:200);
L = length(ttimes);
T_samp = abs(diff(ttimes(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(12*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);

tic
a = d2{1}(:,1:200); a = a(:);
t = repmat(d.xv{1},1,200);
dbz_fine = [];
for j =1:length(a)-200
    [ts, tsi] = sort(t(j:j+199));
    fdata = a(j:j+199)>0.4;
    fdata = fdata(tsi);
    ft = fft(fdata-mean(fdata),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbz_fine(end+1) = ff(mi+1);
end
toc
%%
mdbz = 1e-3*abs(d.scan.data.post_dbz+d.scan.data.pre_dbz)/2;
ttimes = d.xv{1}(1:200);
L = length(ttimes);
T_samp = abs(diff(ttimes(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2); %nyquist frequency
n_alias = floor(mdbz/nyq); % number of times signal is aliased
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);

t = repmat(d.xv{1},1,200);
d3 = reshape(d2{1},200,200,1024);
tnow = tic;
dbz_fine = zeros(size(d3,3),39801);
parfor k = 1:size(d3,3)
a = squeeze(d3(:,:,k)); a = a(:);
tmp = zeros(1,39801);
for j =1:length(a)-200
    [ts, tsi] = sort(t(j:j+199));
    fdata = a(j:j+199)>0.4;
    fdata = fdata(tsi);
    ft = fft(fdata-mean(fdata),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   tmp(j) = ff(mi+1);
end
dbz_fine(k,:) = tmp;
if ~mod(k,100)
    fprintf('done with %i of %i in %.2f seconds \n',k,size(d3,2),toc(tnow));
end
end

%%
tic
%d = ana_avg; 
%d = ana_avg('sm_dBz_phase_raw_140MHz_R_8451');close all
d = ana_avg; close all
d2=anaRawScale(anaRawUnpack(d.scan,d.data),1e-2);
toc % 11s
d3= squeeze(d2{1});
nloop = 256;
d3 = reshape(d3,size(d3,1),nloop,size(d3,2)/nloop);
d3 = reshape(d3,size(d3,1)*size(d3,2),size(d3,3));
dbz_avg = pi*(abs(d.scan.data.pre_dbz)+abs(d.scan.data.post_dbz));
% now (200 x nloop) x nrep
%%
B = 0.5*pi/abs(diff(d.xv{1}(1:2))); %bandwidth
nfreqs = 400;
ffreqs =1e-3*dbz_avg+B./2*linspace(-1,1,nfreqs);
b=.67; a =.25;
nmeas = 200;
ftime = repmat(d.xv{1},1,nloop);
[wmat tmat]=meshgrid(ffreqs,ftime);
bayscoef=a+b*cos(wmat.*tmat);
%figure(1); clf; hold on; c = 'rgbcmyk'; c = [c c c c c c c];
tic;
clear dbzs
for j =1:size(d3,2)
    fdata = 2*(d3(:,j)>0.4)-1;
    fdata2 = 1-repmat(fdata,1,nfreqs).*bayscoef;
    postdist = prod(fdata2(1:nmeas,:)); postdist = postdist/sum(postdist);
    for k = 1:(size(fdata,1)-nmeas)
        postdist = postdist.*fdata2(k+nmeas,:)./fdata2(k,:);
         postdist = postdist/sum(postdist);
         [~,mi]=max(postdist);
         dbzs(j,k) = ffreqs(mi)/2/pi;
         %plot(ffreqs/(2*pi),postdist,c(k)); pause
    end
    if ~mod(j,25)
       fprintf('done with %i of %i in %.2f seconds \n',j,size(d3,2),toc); 
    end
end
toc

%%
ttimes= (1:size(dbzs,2))*4e-6;
dbzs2 = dbzs;
dbzs2(abs(diff(dbzs))>.0001) = NaN;
clear slp
for j = 1:size(dbzs,1)
    mask = ~isnan(dbzs2(j,:));
  tmp = polyfit(ttimes(mask),dbzs2(j,mask),1); 
  %figure(1); clf; hold on
  %plot(ttimes(mask),dbzs2(j,mask));
  %plot(ttimes(mask),polyval(tmp,ttimes(mask)),'r'); pause
  slp(j) = tmp(1);
end

%%
res = [140, -4.87e6; 180, -4.28e6; 30, -4.58e6; 66, -4.85e6; 100, -5.46e6];

%% package stuff for aussies
files =...
    {'sm_dBz_phase_raw_178MHz_R_8452.mat','sm_dBz_phase_raw_140MHz_R_8451.mat',...
    'sm_dBz_phase_raw_30MHz_R_8450.mat', 'sm_dBz_phase_raw_66MHz_R_8449.mat',...
    'sm_dBz_phase_raw_99MHz_R_8448.mat'};


clear results;
tic
for j = 1:length(files)
d = ana_avg(files{j}); close all
d2=anaRawScale(anaRawUnpack(d.scan,d.data),1e-2);
d3= squeeze(d2{1});
nloop = d.scan.data.conf.nloop;
d3 = reshape(d3,size(d3,1),nloop,size(d3,2)/nloop);
results(j).data = d3;
results(j).fname = files{j};
results(j).evotimes = repmat(d.xv{1},1,nloop);
results(j).mdbz = 0.5*(d.scan.data.pre_dbz+d.scan.data.post_dbz);
end
toc

