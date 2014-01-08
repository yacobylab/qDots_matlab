function out = ana_rescale_dbz_v2(file, config)
%function out = ana_rescale_dbz_v2(file, config)
% rescales dbz based on shots taken at a single evolution time to inscrease
% t2*
% This will deal with a true bayesian estimator. broken for now. 

if ~exist('file','var') || isempty(file)
   file = uigetfile('sm*.mat'); 
end

if  ~exist('config','var')
    config = struct();
end


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
T_evo = mode(xv);
%shots = find(xv==T_evo);
%if length(find(diff(shots)>1))==1 % the T_evo also appears in the dbz meas
%    shots(1+find(diff(shots)>1))=[];
%end

config = def(config,'T_rng',[0 Inf]); 

%evo_inds = 1:length(xv); evo_inds(shots) = []; %indicies of evolutions
%use the last consecutive group of evo times 
firstevo = 1+find(abs(diff(xv))>1,1,'last');
evo_inds = firstevo:size(data,2);
evos = xv(evo_inds);  % these are the evolution times
nevo = length(evos);
rrng = 1:nevo; rrng=rrng(rrng>config.T_rng(1)&rrng<config.T_rng(2));
shots = xv; shots(evo_inds) =[];

if isopt(config,'nobayes')
   scale = mean(data(:,shots),2);
else
   scale= findscale(data,xv,shots); 
end
S_av = mean(scale);
%scale = scale-S_av;
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
ti=linspace(1,nevo+100,nevo);
rep=1:length(scale);
alpha=linspace(0.5,2,500);

for i=1:length(alpha);
    for j=1:Nrep;
        %ts=(1+alpha(i)*scale(j))*t;
        ts=((alpha(i))*scale(j))*t;
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
        %keyboard
    end
    dBx(:,i)=mean(dataI);
%     if any(isnan(dBx(:,i)))
%        keyboard 
%     end
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);

%rrng=2:nevo;

fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^config.decay);
c=[];
%first without scaling, to get dbz_avg and deltaDBZ (1/t2*)
beta0 = [.3, .3, 14, .001, 80 0 1.1];
pars = fitwrap(fitopts,evos,mean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 1 0]);
dbz_avg = 2*pi/pars(3);
d_dbz = 1/abs(pars(5));
out.dbzmean=pars; 

%now fit all scales
for j = 1:length(alpha)
    beta0 = pars;%[.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap(fitopts,(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end
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
AMP=.5*exp(-(T_evo/T2)^config.decay);%max(dBx(rrng,I));
guess_alpha = 1/(sqrt((AMP)^2-S_av^2)*T_evo*(dbz_avg+d_dbz));
if ~isreal(guess_alpha)
    fprintf('predicting complex scaling factor. ignoring... \n');
else
    plot(guess_alpha*[1 1], YL,'r');
    guess_alpha = 1/(sqrt((AMP)^2-S_av^2)*T_evo*(dbz_avg-d_dbz));
    plot(guess_alpha*[1 1], YL,'r');
    legend({'Fit Scale Parameter','Calculated Scale Parameter'});
end
title(sprintf('T2*=%.0f',T2));
figind=figind+1;
out.alpha = alpha;
out.t2s = c; 

figure(figind); figind=figind+1; clf; hold on;
keepers = setdiff(1:length(oldscale),ignore);
plot(oldscale,'b.'); plot(keepers,oldscale(keepers),'r.'); 
title(sprintf('signal (mean subtracted), std=%.02f',std(scale)));
xlabel('rep number');
legend({'ignored','included'});

for j=1:Nrep;
    %ts=(1+alphabest*scale(j))*t;
    ts=((alphabest)*scale(j))*t;
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
out.newdata = dBx(:,I);
figure(figind); clf;
%plot(t,sum(data(:,evo_inds))/Nrep,'k.-',talpha(:,I)./(1+alphaV(:,I)),dBx(:,I),'r.-')
plot(t,sum(data(:,evo_inds))/Nrep,'k.-',talpha(:,I),dBx(:,I),'r.-')
%plot(talpha(:,I),dBx(:,I),'r.-')
hold on;

figure(figind); figind = figind+1;
%plot((talpha(rrng,I)./(1+alphaV(rrng,I))),fitfn(pars,(talpha(rrng,I)./(1+alphaV(rrng,I)))));
plot((talpha(rrng,I)),fitfn(pars,(talpha(rrng,I))));
xlabel('evolution time (ns)');
ylabel('Singlet Probability');
title(sprintf('best results: T_{2}^{*} = %.0f',abs(pars(5))));
% hold on
% [c2,gof2] = fit(talpha(rrng,I)./(1+alphaV(rrng,I)),dBx(rrng,I),f,'problem',2);
% plot(c2)
axis([0 nevo 0 0.8]);


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
%bestbetas=bestbetas/pi;

% betas = linspace(-1,1,500);
% probs = zeros(1,length(betas));
% %proball = zeros(size(s.data,1),length(betas));
% MLE = zeros(size(data,1),1);
% bestbeta = MLE;
% for ind = 1:size(data,1)
% for j = 1:length(betas)
%   probs(j) = prod(1+data(ind,:)*betas(j)); 
% end
% probs = probs/sum(probs);
% [m, mi]=max(probs);
% MLE(ind)=m;
% bestbeta(ind) = betas(mi);
% end
% bestbeta= acos(bestbeta);
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

