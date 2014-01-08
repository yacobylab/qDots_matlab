%%
file = uigetfile('*phase*');
s=load(file);
data=s.data;
scan=s.scan;

scantime=getscantime(scan,data);
[t1t t1] = att1('right',scantime,'after'); % FIX ME
data_all=anaHistScale(scan,data,t1); % FIX ME
data=squeeze(data_all{1});
xv=plsinfo('xval',s.scan.data.pulsegroups.name,[],scantime); %FIX ME
xv = round(xv);
firstevo = 1+find(abs(diff(xv))>1,1,'last');
evo_inds = firstevo:size(data,2);
evos = xv(evo_inds); 
% these are the evolution times
nevo = length(evos);
rrng = 1:nevo; %rrng=rrng(rrng>config.T_rng(1)&rrng<config.T_rng(2));
shots = xv; shots(evo_inds) =[];
%%
data2=data;
mm = find(xv==0);
xv(mm)=[];
data2(:,mm)=[];
data2=-1+2*data2;
data2 = 2*(data2>0)-1;
%%
betas = linspace(-.5,.5,2000);
bestbetas = zeros(1,size(data2,1));
%xv = xv/xv(1);
tic;
prob = zeros(1,length(betas));
for j = 1:size(data2,1)
    for b = 1:length(betas)
      prob(b) = prod(1+data2(j,shots).*cos(pi*xv(shots)*betas(b)));
    end
    prob = prob/sum(prob);
    %figure(1); plot(betas,prob); pause
    [m mi]=max(prob);
    bestbetas(j) = betas(mi);
end

%scale = acos(bestbetas)/pi;
toc
scale = (bestbetas);
scale_raw = scale;
%%
S_av = mean(scale);
scale = scale-S_av;
oldscale = scale;
%config = def(config,'s_mask',[min(scale) max(scale)]); % range of points to include
%ignore = find(scale<config.s_mask(1) | scale>config.s_mask(2));
ignore = [];
%data(ignore,:)=[];
%scale(ignore)=[];
%S_av = mean(scale);
%scale = scale-S_av;
out.scale = scale;
out.S_av = S_av;
Nrep = length(scale);

t=linspace(1,nevo,nevo);
ti=linspace(1,nevo+50,nevo);
rep=1:length(scale);
%alpha=1+linspace(-.1,.1,500);
alpha = linspace(-.5,.5,100);
clear tsV; clear dataI;
for i=1:length(alpha);
    for j=1:Nrep;
        ts=(1+alpha(i)*scale(j))*t;
        %ts=(1/scale(j))*t; 
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
        %keyboard
    end
    dBx(:,i)=mean(dataI);
%     if any(isnan(dBx(:,i)))
%        keyboard 
%     end

end
%%
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);

%rrng=2:nevo;

fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
c=[];
%first without scaling, to get dbz_avg and deltaDBZ (1/t2*)
beta0 = [.3, .3, 14, .001, 80 0 1.1];
pars = fitwrap('plinit plfit',evos,mean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 1 0]);
dbz_avg = 2*pi/pars(3);
d_dbz = 1/abs(pars(5));
out.dbzmean=pars; 
%%
%now fit all scales

for j = 1:length(alpha)
    beta0 = pars;%[.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap(fitopts,(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end
%%
fitopts = 'plinit plfit';
c(c>5000)=0;
[T2,I]=max(c);
out.T2 = T2;
alphabest=alpha(I);
figure(figind);  clf;
plot(alpha,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');
%pars = fitwrap('plinit plfit',(talpha(rrng,I)./(1+alphaV(rrng,I)))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);
pars = fitwrap(fitopts,(talpha(rrng,I))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);
figure(figind); hold on;
YL = get(gca,'YLim');
%AMP=.5*exp(-(T_evo/T2)^config.decay);%max(dBx(rrng,I));
%guess_alpha = 1/(sqrt((AMP)^2-S_av^2)*T_evo*(dbz_avg+d_dbz));
% if ~isreal(guess_alpha)
%     fprintf('predicting complex scaling factor. ignoring... \n');
% else
%     plot(guess_alpha*[1 1], YL,'r');
%     guess_alpha = 1/(sqrt((AMP)^2-S_av^2)*T_evo*(dbz_avg-d_dbz));
%     plot(guess_alpha*[1 1], YL,'r');
%     legend({'Fit Scale Parameter','Calculated Scale Parameter'});
% end
title(sprintf('T2*=%.0f',T2));
figind=figind+1;
out.alpha = alpha;
out.t2s = c; 

%%

%% get the files
files = uigetfile('sm*.mat','MultiSelect','on');
if ~iscell(files) 
    files = {files};
end
files = sort(files);

results = [];
ampl=[];
ampfunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));  
eps = zeros(1,length(files));
for j = 1:length(files)
    [a b c]=ana_echo(files{j},struct('opts','echo power rmoutlier','grng', [.21 Inf],'dt',1:64));
    ampl(j,:) = ampfunc(b.params);
    pp=plsinfo('params',c.scan.data.pulsegroups(2).name,[],b.scantime);
    eps(j) = 2*pp(2);
    results(end+1).T2= b.T2;
    results(end).Q = b.Q;
    results(end).T = b.T;
end
%close all;
%%
figure(2); clf;
plot(b.ts(2:end),ampl')
legstr =cellstr(num2str(eps'));
legstr = cellfun(@(x)strcat('eps= ',x),legstr,'UniformOutput',0);
legend(legstr);
xlabel('evolution time (\mus)');
ylabel('echo amplitdue');
%title('two pulse CPMG');

%% get the files
files = uigetfile('sm*.mat','MultiSelect','on');
if ~iscell(files) 
    files = {files};
end
files = sort(files);
legstr = {};
results = [];
ampl=[];
ampfunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));  
eps = zeros(1,length(files));
data = {};
for j = 1:length(files)
    [a b c]=ana_echo(files{j},struct('opts','echo power rmoutlier','grng', [.21 Inf],'dt',1:64));
    ampl(j,:) = ampfunc(b.params);
    pp=plsinfo('params',c.scan.data.pulsegroups(2).name,[],b.scantime);
    eps(j) = 2*pp(2);
    results(end+1).T2= b.T2;
    results(end).Q = b.Q;
    results(end).T = b.T;
    data{j} = squeeze(nanmean(c.data{1}));
    if isempty(str2num(files{j}(14:15)))
        legstr{end+1} = files{j}(14);
    else
       legstr{end+1} = files{j}(14:15); 
    end
    
end
%close all;
figure(2); clf;
plot(b.ts(2:end),ampl')
legend(legstr);
xlabel('evolution time (\mus)');
ylabel('echo amplitdue');

%%
am=[];
for j = 1:length(data)
   am(j,:) =  range(data{j}(2:end,:),2);
end
%%
c='rgbcmykr';c=[c c c c c c c c c];
figure(2); clf; hold on
off = 0;
for j= 1:size(am,1)
  plot(b.ts(2:end),off+am(j,:),c(j))
  off = off +0;.3;%1e-3;
end
legend(legstr);
xlabel('evolution time (\mus)');
ylabel('echo amplitdue');

%% using ana_avg

for j = 1:length(files)
   d=ana_avg(files{j}); 
   data{j} = squeeze(nanmean(d.data{1}));
end