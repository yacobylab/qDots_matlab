%% get files
files = uigetfile('sm*.mat','MultiSelect','on');
if ~iscell(files)
    files = {files};
end
clear out; 

%% analyze single group rot ramsey FID
tic;
for jj = 1:length(files)
figind = 101;
d=ana_avg(files{jj},struct('opts','noppt')); %close all; 
data = squeeze(d.data{1}); 
pp=plsinfo('xval',d.scan.data.pulsegroups.name,[],d.scantime);
xv = pp(2,:);
out(jj).eps = pp(1,1);
out(jj).filename = files{jj};
shots = 1:find(abs(diff(xv))>1e-3);
evo_inds = (shots(end)+1):size(data,2);
scale = squeeze(nanmean(data(:,shots),2));
nanmask = find(isnan(scale));
scale(nanmask) = [];
data(nanmask,:)=[];
evos = xv; evos(shots)=[]; nevo = length(evos);

plot(squeeze(nanmean(d.data{1}(:,:,evo_inds))));
xlabel('evolution time (ns)');
ylabel('P_{singlet}');
title('raw data');


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
Nrep = length(scale);

t=linspace(1,nevo,nevo);
ti=linspace(1,nevo+50,nevo);
rep=1:length(scale);
%alpha=1+linspace(-.1,.1,500);
alpha = linspace(-.6,.6,100);

clear tsV; clear dataI; clear dBx;
for i=1:length(alpha);
    for j=1:Nrep;
        ts=(1+alpha(i)*scale(j))*t;
        %ts=(1/scale(j))*t; 
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
        %keyboard
    end
    dBx(:,i)=nanmean(dataI);
%     if any(isnan(dBx(:,i)))
%        keyboard 
%     end

end

[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);
%rrng=2:nevo;
try
%fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
fitfn = @(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2);
c=[];
%first without scaling, to get dbz_avg and deltaDBZ (1/t2*)
beta0 = [.3, .3, 14, .001, 80 0 1.2];
%pars = fitwrap('',evos,mean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 0 0]);
%pars = fitwrap('plinit plfit',evos,mean(data(:,evo_inds),1),pars,fitfn, [1 1 1 1 1 0]);
%dbz_avg = 2*pi/pars(3);
%d_dbz = 1/abs(pars(5));

pars = fitoscillations(evos,mean(data(:,evo_inds)),'gauss','plinit plfit woff');
dbz_avg = pars(4)/(2*pi);
d_dbz = abs(pars(6));
out(jj).dbzmean=pars; 

%now fit all scales

for j = 1:length(alpha)
    beta0 = pars;%[.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap(fitopts,(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    %pars = fitwrap('',(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 0 0]);
    %pars = fitwrap('plinit plfit',(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]); pause
    %c(j) = abs(pars(5));
    pars=  fitoscillations((talpha(rrng,j))',dBx(rrng,j)','gauss','woff');
    c(j) = abs(1/pars(6));
    
    %pause
end

fitopts = 'plinit plfit';
c(c>1000)=0;
[T2,I]=max(c);
out(jj).T2 = T2;
alphabest=alpha(I);
figure(figind);  clf;
plot(alpha,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');
%pars = fitwrap(fitopts,(talpha(rrng,I))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);

pars = fitoscillations((talpha(rrng,I))',dBx(rrng,I)');
figure(figind); hold on;
YL = get(gca,'YLim');
plot([0 0],YL,'g');
title(sprintf('T2*=%.0f',T2));
figind=figind+1;

out(jj).t2s = c; 

%pars = fitwrap('plinit plfit',(talpha(rrng,I))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);
pars = fitoscillations(talpha(rrng,I)',dBx(rrng,I)','gauss','plinit plfit woff');
figure(figind); clf; hold on;
plot(evos,mean(data(:,evo_inds)),'x-k');
plot((talpha(rrng,I))',dBx(rrng,I)','rx');
plot(talpha(rrng,I),fitfn(pars,(talpha(rrng,I))'),'r');
legend({'unscaled', 'scaled','scaled fit'});
title(sprintf('scaled T2* = %.2f \n unscale T2* = %.2f',T2,abs(1/d_dbz)));
close(500);
out(jj).bareT2 = abs(1/d_dbz);
catch
   out(jj).T2 = NaN;
   out(jj).t2s = NaN(1,200);
   out(jj).dbzmean = NaN(1,6);
   out(jj).bareT2 = NaN;
end
out(jj).alpha = alpha;
out(jj).file = files{jj};

end
toc

%%
close all;
figure(1); hold on;
plot([out.eps],[out.T2]./[out.bareT2],'.');
%set(gca,'YLim',[0 10]);
xlabel('\epsilon_{offset} (mV)');
ylabel('Improvement in T2*');

%%
for j = 1:length(out)
   figure(1); clf; 
   plot(out(j).alpha,out(j).t2s);
   title(sprintf('eps = %.2f',out(j).eps))
   pause
end

%%
d=ana_avg; close all;
data = d.data{1};
pp=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
xv = pp(2,:);
evo_inds = 1+find(abs(diff(xv))>1e-4);
shots = 1:length(xv); shots(evo_inds) = [];
dbzs = squeeze(nanmean(data(:,:,shots),3));

[N X]=hist(dbzs(:),30);
figure(1); clf; hold on;
plot(X,N);
mdbz = nanmean(dbzs(:));
sigma_dbz = nanstd(dbzs(:));
YL = get(gca,'Ylim');
plot(mdbz*[1 1],YL,'g');
plot((mdbz+sigma_dbz*[1 1]),YL,'r')
plot((mdbz+sigma_dbz*[-1 -1]),YL,'r')


thold = 10;
mask = ~(dbzs>mdbz-thold*sigma_dbz & dbzs < mdbz+thold*sigma_dbz);
datatmp = reshape(data(:,:,evo_inds),size(data,1)*size(data,2),length(evo_inds));
datatmp(mask(:),:) = NaN;
datatmp = reshape(datatmp,size(data,1),size(data,2),length(evo_inds));
datatmp = squeeze(nanmean(datatmp));
figure(2); clf;
imagesc(datatmp);
dbztmp = dbzs(:);
figure(3); clf; hold on;
plot(dbzs(:),'.');
dbztmp(mask(:))=[];
plot(dbztmp,'r.');
legend({'excluded','kept'});
title(sprintf('keeping %.1f%% of the dbz values',100*numel(dbztmp)/numel(dbzs)));


%%
Ts= linspace(.3,3,20); %HACK!
ampls=[]; T2a = [];
avals = [10 5 2 1 .8 .5 .4 .3 .2 .1 .05 .01];
for j = 1:length(avals)
thold = avals(j);
mask = ~(dbzs>mdbz-thold*sigma_dbz & dbzs < mdbz+thold*sigma_dbz);
datatmp = reshape(data(:,:,evo_inds),size(data,1)*size(data,2),length(evo_inds));
datatmp(mask(:),:) = NaN;
datatmp = reshape(datatmp,size(data,1),size(data,2),length(evo_inds));
datatmp = squeeze(nanmean(datatmp));
figure(2); clf;
imagesc(datatmp);
dbztmp = dbzs(:);
figure(3); clf; hold on;
plot(dbzs(:),'.');
dbztmp(mask(:))=[];
plot(dbztmp,'r.');
legend({'excluded','kept'});
title(sprintf('keeping %.1f%% of the dbz values',100*numel(dbztmp)/numel(dbzs)));
ampls(:,j) = range(datatmp,2);
P=polyfit(Ts',ampls(:,j),1);
T2a(j) = 1/abs(P(1));
pause
end

c = 'rgbcmyk'; c= [c c c  c];
figure(4); clf; hold on;
plot(ampls)
legstr =cellstr(num2str(avals'));
legstr = cellfun(@(x)strcat('threshold = ',x),legstr,'UniformOutput',0);
legend(legstr);
