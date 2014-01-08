function out = ana_rescale_dbz_v3(file, config)
%function out = ana_rescale_dbz_v2(file, config)
% rescales dbz based on shots taken at a single evolution time to inscrease
% t2*
% This will deal with a true bayesian estimator. broken for now. 

if ~exist('file','var') || isempty(file)
   file = uigetfile('sm*.mat'); 
end
if isempty(file)
    return
end

if  ~exist('config','var')
    config = struct();
end

out.file = file;
config = def(config,'opts','');
config = def(config,'decay',1.5);
config = def(config,'rng',[1 inf]);
config = def(config,'side','right');
config = def(config,'figind',470); figind = config.figind; fb = figind;

s=load(file);
data=s.data;
scan=s.scan;

scantime=getscantime(scan,data);
[t1t t1] = att1(config.side,scantime,'after'); % FIX ME
data_all=anaHistScale(scan,data,t1); % FIX ME
data=squeeze(data_all{1});

if isopt(config,'plotfit')
    fitopts = 'plinit plfit';
else
    fitopts = '';
end
config = def(config,'frames',1:size(data,1)); frames = config.frames;
data = data(frames,:);

% determine which shots are used for scaling.
xv=plsinfo('xval',s.scan.data.pulsegroups.name,[],scantime); %FIX ME
if size(xv,1)>1
    xv = xv(2,:);
end
config = def(config,'evo_inds',find(diff(xv)==1)); evo_inds = [config.evo_inds, config.evo_inds(end)+1];
%evo_inds = find(diff(xv)==1); evo_inds = [evo_inds, evo_inds(end)+1];
evos = xv(evo_inds);
shot_inds = 1:length(xv);
shot_inds(evo_inds) = [];
shots = xv(shot_inds);

[shot_sort, inds]=sort(shots);

num_evos = length(unique(shots)); % bit of a hack
shot_sort = reshape(shot_sort,numel(shot_sort)/num_evos,num_evos);
inds = reshape(inds,numel(inds)/num_evos,num_evos);

if 1 %avg together shots of same evotime
%shot_data = data(:,shot_inds);
%shot_data = 1-2*squeeze(mean(reshape(shot_data,size(shot_data,1),size(shot_sort,1),size(shot_sort,2)),2));
evotime = shots(inds(1,:));
for j = 1:size(data,1)
   shot_data(j,:) = mean(reshape(data(j,shot_inds),num_evos,(length(shot_inds)/num_evos)),2); 
   shot_data(j,:) = 1-2*shot_data(j,:)>0.5;
end
else
shot_data = 1-2*data(:,shot_inds)>0;
%shot_data = 1-2*squeeze(mean(reshape(shot_data,size(shot_data,1),size(shot_sort,1),size(shot_sort,2)),2));
evotime = shots(inds(:));
end
config = def(config,'shot_include',1:size(shot_data,2));
shot_data = shot_data(:,config.shot_include);
%now shot_data is now #repitions x # different evo times per line
% evo time is the different evo times


config = def(config,'T_rng',[0 Inf]); 

%first get the raw data and analyze
rawdata = mean(data(:,evo_inds));
rawpars = fitdbzoscillations(evos,rawdata,'plinit,plfit');
baret2 = 1/rawpars(6); %in ns
omega_dbz = rawpars(4)/(2*pi); %in GHz
config = def(config,'mdbz',omega_dbz'); % the mean dbz
out.f_dbz = omega_dbz;
fprintf('dbz = %d MHz \n bare T2* = %d ns \n', omega_dbz*1e3,baret2);
out.rawpars = rawpars;

nevo = length(evo_inds);
%t=linspace(1,nevo,nevo);
%ti=linspace(1,nevo+50,nevo);
t=linspace(evos(1),evos(end),nevo);
ti=linspace(evos(1),evos(end)+50,nevo);

Nrep = size(data,1); rep = 1:Nrep;
if 1;%numel(evotime)==1
omega=(omega_dbz+3*(1/baret2)*linspace(-1,1,500));
else
%omega=(omega_dbz+0.9*abs(1/(evotime(2)-evotime(1)))*linspace(-1,1,500));
%dd = diff(evotime);
%dd = dd(find(dd>0,1));
%omega = omega_dbz + (.5/dd)*linspace(-1,1,500);
omega = linspace(0,omega_dbz+5/baret2,1000);
end
%see matlab documentation for FFT normalization
L = length(evotime);
T_samp = abs(diff(evotime(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(6*L);
freqs = getfreqs(F_samp, NFFT,config.mdbz);
%freqs = F_samp-(F_samp/2)*linspace(0,1,NFFT/2+1);%omega_dbz+(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
prior = exp(-(freqs-omega_dbz).^2/(2*rawpars(6)^2));prior = prior/sum(prior);
%ts = (0:L-1)*T_samp;
for j = 1:Nrep
    if 0
    Ndiv = 4; np = size(shot_data,2)/Ndiv;
    ff=0;
    for ii= 1:Ndiv
        NFT2 = 2^nextpow2(10*np);
        mask =(1+(ii-1)*np):ii*np;
        ft = fft(shot_data(j,mask)-mean(shot_data(j,mask)),NFT2);
        ff = ff+abs(ft(1:NFT2/2+1));
    end
    fr2=linspace(freqs(1),freqs(end),NFT2/2+1);
    [~, mi]=max(ff);
    pmax(j) = fr2(mi);
    %figure(1); clf; plot(fr2,ff); keyboard
    else
        if isopt(config,'flip')
           mask = size(shot_data,2):-1:1; 
        else
            mask = 1:1:size(shot_data,2);
        end
    ft = fft(shot_data(j,mask)-mean(shot_data(j,mask)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    [~, mi] = max(ft(2:end));
    pmax(j) = freqs(mi+1); %+1 bc 2:end
    end
    if 0
    f2 = freqs;
    for k = 1:length(f2)
       %p2max(k) = abs(trapz(shot_data(j,:).*exp(1i*2*pi*f2(k)*evotime)));
       p2max(k) = (trapz(shot_data(j,:).*cos(2*pi*f2(k)*evotime)));
    end
    p2max(1:4)=0; p2max(end-3:end)=0;
    %figure(1); clf; plotyy(freqs,ft,f2,p2max); pause
    %keyboard
    [~,mi]=max(p2max(2:end-1));
     pmax(j) = f2(mi);
    end
    if 0
        figure(1); clf;  hold on
        plot(freqs,ft);
        YL = get(gca,'YLim');
        plot(pmax(j)*[1 1],YL,'r');
        keyboard
    end
    if 0
        omega = linspace(min(freqs),max(freqs),length(freqs));
    %p_omega = exp(-(omega-omega_dbz).^2*(2*(baret2*pi)^2)); p_omega(1:10)=0;
    p_omega = ones(size(omega));
    p_omega = p_omega/sum(p_omega);
    %figure(1); plot(omega,p_omega); pause
    %p_omega = ones(1,length(omega));
    V = 1;%range(rawdata);
    for k = 1:size(shot_data,2)
       p_omega = p_omega.*(1+shot_data(j,k)*cos(evotime(k)*omega));
       %p_omega = p_omega.*(1+V*exp(-evotime(k)^2/baret2^2)*shot_data(j,k)*cos(evotime(k)*omega));
       %figure(1); plot(omega,p_omega);keyboard;
    end
    %p_omega(1:100) = -Inf; p_omega(end) = -Inf;
   [m mi]=max(p_omega);
   p2max(j) = omega(mi);
  %figure(1); plot(omega,p_omega);pause;
   %keyboard
    end
    if 0
       omega_fine = pmax(j)+ .4*diff(freqs(1:2))*linspace(-.5,.5,100);
       p_omega = ones(size(omega_fine));p_omega = p_omega/sum(p_omega);
       for k = 1:size(shot_data,2)
          p_omega = p_omega.*(1+shot_data(j,k)*cos(evotime(k)*omega_fine))/sum(p_omega);
          %figure(1); clf; plot(omega_fine,p_omega); pause
       end
       %figure(1); clf; plot(omega_fine,p_omega); hold on; YL = get(gca,'YLim'); plot(pmax(j)*[1 1],YL,'r');pause
       [~,mi] = max(p_omega);
       if mi >0 && mi < length(p_omega)
        pmax(j) = omega_fine(mi);
       end
    end
end
if exist('p2max','var')
    out.p2max = p2max;
end
%tfunc = @(t,a,s) t/(1+a*s); alpha = linspace(-1,1,500);
%tfunc = @(t,a,s) t*(1/(a*(s+.5))); alpha = fudge*linspace(.3,1.7,500);%1/fudge+linspace(-.1,.1,500);
%tfunc = @(t,a,s) t*(1/((a+s+.5))); alpha = .5*linspace(-1,1,500);%1/fudge+linspace(-.1,.1,500);
tfunc = @(t,a,s) t*(((s+.5))); alpha = .5*linspace(-.5,.5,500);%1/fudge+linspace(-.1,.1,500);
out.tfunc = tfunc;

figure(figind); figind = figind+1;
plot(pmax,'.'); title('all values of estimated dbz');
[h b]=hist(pmax,40);
figure(figind); figind = figind+1;
plot(b,h); 
xlabel('\Delta{}B_z'); ylabel('counts');
scale = -.5+pmax/omega_dbz;
out.scale = scale;
out.pmax = pmax;
fudge = omega_dbz/mean(pmax);
out.fudge = fudge;
out.data = data;
fitfn = @(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2);
for j = 1:Nrep
   %datatmp(j,:) = qinterp1(t*fudge*(pmax(j)/omega_dbz)^-1,data(j,evo_inds),ti);
   datatmp(j,:) = qinterp1(tfunc(t,0,scale(j)),data(j,evo_inds),ti);
end
datascaled = nanmean(datatmp);
mask = ~isnan(datascaled);
% if sum(mask)==0
%    error('interpolated data is all nans'); 
% end
if isopt(config,'quick')
    return
end
%try
scaled_pars = fitdbzoscillations(evos(mask),datascaled(mask),'plinit plfit');
figure(figind); clf; figind = figind+1;
hold on;
plot(evos,rawdata,'kx-');
plot(evos,datascaled,'r.');
plot(evos,fitfn(scaled_pars,evos));
xlabel('evolution time (ns)'); ylabel('Singlet Probability');
title(sprintf('best results no fitting: T_{2}^{*} = %.0f',1/abs(scaled_pars(6))));
%catch
%   fprintf('cannot fit scaled data'); 
%end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);
rrng = 1:nevo; rrng=rrng(rrng>config.T_rng(1)&rrng<config.T_rng(2));
if ~isopt(config,'fitscale') %change this to actually fit the scale
    alphabest = fudge;
    fitscale = 0;
else
    fitscale = 1;
for i=1:length(alpha);
    for j=1:Nrep;
        ts = tfunc(t,alpha(i),scale(j));
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
        %keyboard
    end
    dBx(:,i)=nanmean(dataI);
end
%now fit all scales
for j = 1:length(alpha)
   %figure(33); plot(talpha(rrng,j),dBx(rrng,j)); pause(.005);
   mask = ~isnan(dBx(rrng,j));
    pars = fitdbzoscillations((talpha(mask,j))',dBx(mask,j)','plinit plfit');
    c(j) = 1/abs(pars(6));
    Q(j) = abs(pars(4))/(2*pi*abs(pars(6)));
    %pause
end
c(c>1000)=0;
Q(c>1000)=0;
Q(Q>100)=0;
figure(figind);  clf; figind = figind+1;
plotyy(alpha,c,alpha,Q);
[T2,I]=max(Q);
out.T2 = T2;
alphabest=alpha(I);
out.alphabest = alphabest;
fprintf('best scale factor = %f\n',alphabest);
figure(figind);  clf; figind = figind+1;
plot(alpha,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');
%pars = fitwrap('plinit plfit',(talpha(rrng,I)./(1+alphaV(rrng,I)))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);
%pars = fitwrap(fitopts,(talpha(rrng,I))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);
%mask = ~isnan(dBx(rrng,I));
%pars = fitdbzoscillations((talpha(mask,I))',dBx(mask,I)','plinit plfit');
%figure(figind); hold on;
%figind=figind+1;
out.alpha = alpha;
out.t2s = c; 
out.Qs = Q;
out.newdata = dBx(:,I);
end %

for j=1:Nrep;
    ts = tfunc(t,alphabest,scale(j));
    %ts=(1+alphabest*scale(j))*t;
    %ts=t*(omega_dbz/(alphabest*pmax(j)));
    %ts=t/(1+alphabest*scale(j));
    tsV(j,:)=ts;
    dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti); %in the path of z:/qDots/matlab
end


figure(figind); figind=figind+1;
pcolor(tV,repV,data(:,evo_inds));shading flat;caxis([-1 2])
title('without scaling');
xlabel('evolution time (ns)');
ylabel('reps');
figure(figind); figind=figind+1; clf;
pcolor(tiV,repV,dataI);shading flat; %caxis([-1 2])
title('with scaling');
xlabel('evolution time (ns)');
ylabel('reps');

out.rawdata = mean(data(:,evo_inds));

if fitscale
mask = ~isnan(dBx(rrng,I));
pars = fitdbzoscillations((talpha(mask,I))',dBx(mask,I)','plinit plfit');
figure(figind); clf; hold on;
plot(t,mean(data(:,evo_inds)),'k.-',talpha(:,I),dBx(:,I),'r.-')
plot(t,fitfn(pars,t));
else
   mask = ~isnan(datascaled);
   pars = fitdbzoscillations(t(mask),datascaled(mask),'plinit plfit'); 
   figure(figind); clf; hold on;
   plot(t,out.rawdata,'k.-',t,datascaled,'r.-');
   plot(t,fitfn(pars,t));
end
figind = figind+1;
xlabel('evolution time (ns)');
ylabel('Singlet Probability');
title(sprintf('best results: T_{2}^{*} = %.0f',1/abs(pars(6))));
axis([0 nevo 0 0.8]);
fprintf('Raw Q = %f \n Scaled Q = %f',abs(baret2*omega_dbz),abs(pars(4)/(pars(6)*2*pi)));


if 0
%now plot the best ones

%figure(figind); clf;
%plot(t,sum(data(:,evo_inds))/Nrep,'k.-',talpha(:,I)./(1+alphaV(:,I)),dBx(:,I),'r.-')
%plot(t,sum(data(:,evo_inds))/Nrep,'k.-',talpha(:,I),dBx(:,I),'r.-')
%plot(talpha(:,I),dBx(:,I),'r.-')
%plot(t,out.rawdata,'k.-',t,datascaled,'r.-');
%hold on;

%figure(figind); figind = figind+1;

%plot((talpha(rrng,I)./(1+alphaV(rrng,I))),fitfn(pars,(talpha(rrng,I)./(1+alphaV(rrng,I)))));

%plot(t,fitfn(pars,t));
%plot((talpha(rrng,I)),fitfn(pars,(talpha(rrng,I))));
end
if ~isopt(config,'noppt')
    prettyfile=regexprep(file,'(sm_)|(\.mat)','');
    ppt=guidata(pptplot);
    set(ppt.e_file,'String',file);
    set(ppt.e_figures,'String',['[',sprintf('%d ',fb:figind-1),']']);
    set(ppt.e_title,'String',prettyfile);
    set(ppt.exported,'Value',0);
end
end

function bestbetas = findscale(data,xv,mask)
%renormalize and digitize
mm = find(xv==0);
xv(mm)=[];
data(:,mm)=[];
data=-1+2*data;
data = 2*(data>0)-1;
betas = pi*linspace(.7,1.3,2000);
bestbetas = zeros(1,size(data,1));
xv = xv/xv(1);
for j = 1:size(data,1)
    for b = 1:length(betas)
    prob(b) = prod(1+data(j,mask).*cos(xv(mask)*betas(b)));
    end
    prob = prob/sum(prob);
    [m mi]=max(prob);
    bestbetas(j) = betas(mi);
end
bestbetas = acos(bestbetas)/pi;

 end
% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
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

%evo_inds = 1:length(xv); evo_inds(shots) = []; %indicies of evolutions
%use the last consecutive group of evo times 
% firstevo = 1+find(abs(diff(xv))>1,1,'last');
% evo_inds = firstevo:size(data,2);
% evos = xv(evo_inds);  % these are the evolution times
% nevo = length(evos);
% rrng = 1:nevo; rrng=rrng(rrng>config.T_rng(1)&rrng<config.T_rng(2));
% shots = xv; shots(evo_inds) =[];

% if isopt(config,'nobayes')
%    scale = mean(data(:,shots),2);
% else
%    scale= findscale(data,xv,shots); 
% end
% S_av = mean(scale);
%scale = scale-S_av;
%oldscale = scale;
%config = def(config,'s_mask',[min(scale) max(scale)]); % range of points to include
%ignore = find(scale<config.s_mask(1) | scale>config.s_mask(2));
%ignore = [];
%data(ignore,:)=[];
%scale(ignore)=[];
%S_av = mean(scale);
%scale = scale-S_av;
%out.scale = scale;
%out.S_av = S_av;
%Nrep = length(scale);
