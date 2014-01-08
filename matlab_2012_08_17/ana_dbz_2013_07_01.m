%%
files = get_files('sm*.mat');

%%
results=[];
for i = 1:length(files)
d = ana_avg(files{i});

data = squeeze(d.data{1});
xv = d.xv{1};
evo_inds = find(diff(xv)==1); evo_inds = [evo_inds, evo_inds(end)+1)];
for j = 1:size(data,1)
   shot_data(j,:) = mean(reshape(data(j,shot_inds),num_evos,(length(shot_inds)/num_evos)),2); 
end
evos = xv(evo_inds);
shot_inds = 1:length(xv);
shot_inds(evo_inds) = [];
shots = xv(shot_inds);
num_evos = length(unique(shots)); % bit of a hack
shot_sort = reshape(shot_sort,numel(shot_sort)/num_evos,num_evos);
inds = reshape(inds,numel(inds)/num_evos,num_evos);
Nrep = size(d.data{1},1); rep = 1:Nrep;
evotime = shots(inds(1,:));
L = length(evotime);
T_samp = abs(diff(evotime(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(L);
freqs = omega_dbz+(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
prior = exp(-(freqs-omega_dbz).^2/(2*rawpars(6)^2));prior = prior/sum(prior);
%ts = (0:L-1)*T_samp;
for j = 1:Nrep
    ft = fft(shot_data(j,:)-mean(shot_data(j,:)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    [~, mi] = max(ft(2:end));
    pmax(j) = freqs(mi+1); %+1 bc 2:end

end
scale = -.5+pmax/omega_dbz;
fudge = omega_dbz/mean(pmax);


end

%% a seed one
d_good = ana_rescale_dbz_v3([],struct('opts','fitscale'));
[h b]= hist(d_good.scale,50);
fitfn = @(p,x) p(1)*exp(-((x-p(2))/(p(3))).^2);
pars = fitwrap('plinit plfit',b,h,[max(h), b(end/2) .025],fitfn);
m_good = pars(2);

%% analyze rest
clear results dataI
alphabest = 0;%d_good.alphabest;
four_inds = 1:200;
evo_inds = (1+four_inds(end)):250;
fitfn = @(p,x) p(1)*exp(-((x-p(2))/(p(3))).^2);
for j = 1:length(files)
    d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt flip','mdbz',.14,'shot_include',1:200));
    %d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt ','mdbz',.089));
    results(j).tfunc = d2.tfunc;
    results(j).dbz = d2.f_dbz;
    if 0
        [h b]=hist(d2.scale,50);
        [~, mi]=max(h); 
        pars = fitwrap('plinit plfit',b,h,[max(h), b(mi) .025],fitfn);
        results(j).scale = d2.scale*(pars(2)/m_good);
    else
        %results(j).scale = d2.scale*(d_good.f_dbz/d2.f_dbz);
        results(j).scale = d2.scale;%*(d_good.f_dbz/d2.f_dbz);
    end
    results(j).fudge = d2.fudge;
    results(j).pmax = d2.pmax;
    d=ana_avg(files{j},'noppt'); 
    results(j).file = d.filename;
    results(j).data = squeeze(d.data{1});
    data = squeeze(d.data{1});
    t = d.xv{1}(evo_inds);
    ti = linspace(.95*t(1),t(end)*1.05,length(t));
    %ph_good = linspace(min(t(1)*results(j).pmax)-pi,pi+max(t(end)*results(j).pmax),length(t));
    %factor = 5;
    for k = 1:size(data,1)        
        %tnew = linspace(ts(1),ts(end),length(ts)*factor);
        %ph = results(j).scale(k)*t;
        %dataI(k,:)=interp1(ph,data(k,evo_inds),ph_good); 
        %ts = results(j).tfunc(t,0,results(j).scale(k));
        ts = t*results(j).pmax(k)/mean(results(j).pmax);
        if all(ts==0)
           ts = ti; 
        end
        dataI(k,:)=interp1(ts,data(k,evo_inds),ti); 
    end
    mmask = 2:size(dataI,2);
    results(j).newdata = nanmean(dataI(:,mmask));
    results(j).xv = ti(mmask);
    %keyboard
    fprintf('done with %i of %i \n',j,length(files));
end
close all;
%% make plot
figure(3); clf; hold on;
for j = 1:length(results)
    mask2 = 2:length(results(1).xv);
    plot(results(j).xv(mask2),results(j).newdata(mask2));    
end
xlabel('Evolution Time (ns)'); ylabel('Triplet Probability');
title(sprintf('dt = %d',abs(diff(d.xv{1}(1:2)))));
%%
data_all = [results.newdata]; etime = [results.xv];
[ee ii]=sort(etime);
data_all = data_all(ii);
mask = (data_all)<.7 & data_all > -.01;
beta0 = [.36, .34, 12, .01, 1900 0];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^2);
pars = fitwrap('plinit plfit',ee(mask),data_all(mask),beta0,fitfn, [1 1 1 1 1 0]);

figure(4); clf; hold on;
plot(ee(mask),data_all(mask),'k-x'); axis([0 1])
%plot(ee(mask),fitfn(pars,ee(mask)),'r');


%% set up mfit to fit t2* and ignore phase

data = [];
model = [];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
mask = [1 1 1 1 0 0];
beta0 = [.36 .36 12 1400 0 1.9];
offset = 0;
for j = 1:length(results)
    if all((~isnan(results(j).newdata)))
        data(end+1).y = results(j).newdata;
        data(end).x = results(j).xv;
        model(end+1).fn = fitfn;
        model(end).pt = @(p) [p(1:6), p(7+offset)];
        offset = offset+1;
        mask = [mask 1];
        beta0 = [beta0, .01]; 
    end
end

pars = mfitwrap(data,model,beta0,'plinit plfit',mask);
mask(6) = 1;
pars = mfitwrap(data,model,pars,'plinit plfit',mask);



%% plot the mfit
figure(44); clf; hold on;
for j = 1:length(data)
    plot(data(j).x,data(j).y,'x');
    plot(data(j).x,fitfn(model(j).pt(pars),data(j).x),'k');
end
xlabel('Evolution Time (ns)');
ylabel('P_{Triplet}');
title(sprintf('T_2^* = %.1f ns \n decay exponent = %.1f',pars([4,6])));
%% a nice string for the ppt
pptstr = sprintf('dt = %d \n',abs(diff(d.xv{1}(1:2))));
pptstr = [pptstr sprintf('%s ',[files{1} ' to ' files{end}])];
clipboard('copy', pptstr);
pptplot

%% 
files = get_files('sm*.mat');
%% analyze the sideways data
evo_inds = 201:250;
four_inds = 1:201;
mdbz = .089;
all_phase = [];
all_data = [];
all_time = [];
errs = [];
tic
for j = 1:length(files)
d= ana_avg(files{j},struct('opts','noplot noppt'));
four_data = squeeze(d.data{1}(:,:,four_inds));
four_data = four_data>0.5;
evo_data = mean(squeeze(d.data{1}(:,:,evo_inds)),2);
%evo_std = std(squeeze(d.data{1}(:,:,evo_inds)),[],2)';errs = [errs,evo_std];
L = size(four_data,2);
NFFT = 2^nextpow2(6*L);
xv = plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
Tsamp = abs(xv(2,1)-xv(2,2)); Fsamp=1/Tsamp;
%Tsamp = diff(abs(d.xv{1}(four_inds(1:2)))); Fsamp = 1/Tsamp;
nyq = abs(Fsamp/2);
n_alias = floor(mdbz/nyq);
ff = (n_alias+.5)*nyq+(-1)^n_alias*(Fsamp/2)*linspace(-.5,.5,NFFT/2+1);
ph=[]; frq=[]; etimes = [];
for k = 1:size(four_data,2)
    for ll = 1:size(four_data,1)
    ft = fft(four_data(ll,k,:)-mean(four_data(ll,k,:)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));
    [~, mi] = max(ft(2:end));
    %figure(1); clf; plot(ff,ft); pause
    frq(ll) = ff(mi+1); %+1 bc 2:end
    %etimes(k) = d.xv{1}(evo_inds(1))-(28-d.xv{k}(evo_inds(1)));
    etimes(ll) = d.xv{k}(evo_inds(1));
    ph = 2*pi*frq.*etimes;
    all_phase = [all_phase, ph];
    all_data = [all_data, squeeze(evo_data)'];
    all_time = [all_time, etimes];
    end
end
%etimes(1) = etimes(2)-1;
j



end
toc
%%
[all_phase2,inds]=sort(all_phase);
all_data2=all_data(inds);


%% analyze sweep of how many pts used for estimation

evo_list = [10 25 50 100 125 150 175 200];
clear results_all;
for ll = 1:length(evo_list)
clear results dataI
four_inds = 1:200;
evo_inds = (1+four_inds(end)):250;
for j = 1:length(files)
    d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt flip','mdbz',.089,'shot_include',1:evo_list(ll)));
    %d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt ','mdbz',.089));
    results(j).dbz = d2.f_dbz;
    results(j).scale = d2.scale;
    results(j).pmax = d2.pmax;
    d=ana_avg(files{j},'noppt'); 
    results(j).file = d.filename;
    results(j).data = squeeze(d.data{1});
    data = squeeze(d.data{1});
    t = d.xv{1}(evo_inds);
    ti = linspace(.95*t(1),t(end)*1.05,length(t));
    %ph_good = linspace(min(t(1)*results(j).pmax)-pi,pi+max(t(end)*results(j).pmax),length(t));
    %factor = 5;
    for k = 1:size(data,1)        
        ts = t*results(j).pmax(k)/mean(results(j).pmax);
        if all(ts==0)
           ts = ti; 
        end
        dataI(k,:)=interp1(ts,data(k,evo_inds),ti); 
    end
    mmask = 2:size(dataI,2);
    results(j).newdata = nanmean(dataI(:,mmask));
    results(j).xv = ti(mmask);
    %keyboard
    fprintf('done with %i of %i \n',j,length(files));
end

data = [];
model = [];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
mask = [1 1 1 1 0 0];
beta0 = [.36 .36 12 1400 0 1.9];
offset = 0;
for j = 1:length(results)
    if all((~isnan(results(j).newdata)))
        data(end+1).y = results(j).newdata;
        data(end).x = results(j).xv;
        model(end+1).fn = fitfn;
        model(end).pt = @(p) [p(1:6), p(7+offset)];
        offset = offset+1;
        mask = [mask 1];
        beta0 = [beta0, .01]; 
    end
end
try
pars = mfitwrap(data,model,beta0,'plinit plfit',mask);
mask(6) = 1;
pars = mfitwrap(data,model,pars,'plinit plfit',mask);
catch
pars = beta0;
end
results_all(ll).data = data;
results_all(ll).model = model;
results_all(ll).fitfn = fitfn;
results_all(ll).fitpars = pars;
results_all(ll).t2=pars(4);
results_all(ll).n_evos = evo_list(ll);
end

%% plot stuff

figure(1); clf; 
plot([results_all.n_evos],[results_all.t2],'x','LineWidth',3,'MarkerSize',7)
xlabel('Number of Evolution Times');
ylabel('T_2^* (ns)');
title([sprintf('dt = %d ns\n', abs(diff(d.xv{1}(1:2)))), ' 8 \mu{}s pulse'])

%%
finds_1 = 1:25;
finds_2 = 26:100;
evo_inds = 101:150;
mdbz = .089;
results = [];
for j = 1:length(files)
   d=ana_avg(files{j},'noppt');
   data = squeeze(d.data{1});
   Tsamp = abs(diff(d.xv{1}(finds_1(1:2)))); Fsamp = 1/Tsamp;
   L = length(finds_1);
   NFFT = 2^nextpow2(6*L);
   nyq = abs(Fsamp/2);
   n_alias = floor(mdbz/nyq);
   ff = (n_alias+.5)*nyq+(-1)^n_alias*(Fsamp/2)*linspace(-.5,.5,NFFT/2+1);
   for k= 1:size(data,1)
    ft = fft(data(k,finds_1)-mean(data(k,finds_1)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    [~, mi] = max(ft(2:end));
    pmax(k) = ff(mi+1); %+1 bc 2:end
    %keyboard
   end
   results(j).pmax = pmax;
   if 1
   for k = 1:size(data,1)
    Tsamp = abs(diff(d.xv{1}(finds_2(1:2)))); Fsamp = 1/Tsamp;
    L = length(finds_2);
    NFFT = 2^nextpow2(6*L);
    nyq = abs(Fsamp/2);
    n_alias = floor(pmax(k)/nyq);
    ff = (n_alias+.5)*nyq+(-1)^n_alias*(Fsamp/2)*linspace(-.5,.5,NFFT/2+1);
    ft = fft(data(k,finds_2)-mean(data(k,finds_2)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    [~, mi] = max(ft(2:end));
    p2max(k) = ff(mi+1); %+1 bc 2:end
    %keyboard
   end
   else
      p2max = pmax; 
   end
   results(j).p2max = p2max;
   
   t = d.xv{1}(evo_inds);
   ti = linspace(.95*t(1),t(end)*1.05,length(t));
    for k = 1:size(data,1)        
        ts = t*results(j).p2max(k)/mean(results(j).p2max);
        if all(ts==0)
           ts = ti; 
        end
        dataI(k,:)=interp1(ts,data(k,evo_inds),ti); 
    end
    results(j).newdata = nanmean(dataI);
    results(j).tvals = ti;
end


%% make plot
figure(33); clf; hold on;
for j = 1:length(results)
    mask2 = 2:length(results(1).tvals);
    plot(results(j).tvals(mask2),results(j).newdata(mask2));    
end
xlabel('Evolution Time (ns)'); ylabel('Triplet Probability');
title(sprintf('dt = %d',abs(diff(d.xv{1}(1:2)))));