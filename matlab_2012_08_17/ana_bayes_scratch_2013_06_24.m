w0 = .082*2*pi;
t = 2:2:2*100;
s = cos(w0*t);

omegas = linspace(-.01,.01,100)+w0;
for j = 1:length(omegas)
   C(j) =  trapz(s.*cos(t*omegas(j)));
end

figure(1); clf; plot(omegas/(2*pi),C);

ft = fftshift(fft(s)); ft = abs(ft);
figure(2); clf; 
freqs = linspace(0,1/4,100);
plot(freqs,ft)

%%
d=ana_avg; close all
data = squeeze(d.data{1});
xv = d.xv{1};

% shot_inds = 1:400;
% four_inds = 401:1000;
% evo_inds = 1001:length(xv);

shot_inds = 801:1000;
four_inds = 801:1200;
evo_inds = 1201:length(xv);

shot_data = data(:,shot_inds);
four_data = data(:,four_inds);
evo_data = data(:,evo_inds);

figure(1); clf; 
plot(xv(evo_inds),mean(evo_data));


%% fit raw data
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^2);
c=[];
%first without scaling, to get dbz_avg and deltaDBZ (1/t2*)
beta0 = [.3, .3, 10, .001, 80 0];
pars = fitwrap('plinit plfit',xv(evo_inds),mean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 1 0]);
dbz_avg = 2*pi/pars(3);
d_dbz = 1/abs(pars(5));


%% use the shots at one evo time, first calculate the scale factors
data_dig=-1+2*shot_data;
data_dig = 2*(data_dig>0)-1;
betas = linspace(-1,1,500);
probs = zeros(1,length(betas));
%proball = zeros(size(s.data,1),length(betas));
MLE = zeros(size(data_dig,1),1);
bestbeta = MLE;
for ind = 1:size(data_dig,1)
for j = 1:length(betas)
  probs(j) = prod(1+data_dig(ind,:)*betas(j)); 
end
probs = probs/sum(probs);
[m, mi]=max(probs);
MLE(ind)=m;
bestbeta(ind) = betas(mi);
end
scale= acos(bestbeta);
S_av = mean(scale);
scale = scale-S_av;

%% do the interpolation
nevo= length(evo_inds);
Nrep = length(scale);

t=linspace(1,nevo,nevo);
ti=linspace(1,nevo+50,nevo);
rep=1:length(scale);
alpha=linspace(-0.25,.25,200);

for i=1:length(alpha);
    for j=1:Nrep;
        %ts=t/(1+alpha(i)*scale(j));
        %ts=(1+alpha(i)*scale(j))*t;
        ts = tfunc(t,alphabest,newscale(j))*(1+alpha(i)*scale(j));
        
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=mean(dataI);
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);

%rrng=2:nevo;

%% fit scaled data
rrng = 2:size(talpha,1);
for j = 1:length(alpha)
    beta0 = pars;%[.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end

c(c>5000)=0;
[T2,I]=max(c);
out.T2 = T2;
alphabest=alpha(I);
figure(1);  clf;
plot(alpha,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');
pars = fitwrap('plinit plfit',(talpha(rrng,I))',dBx(rrng,I)',beta0,fitfn, [1 1 1 1 1 0]);

omegas=dbz_avg./(1+alphabest*scale);
new_ddbz = 1/T2;
%%
[h b]=hist(omegas,40);
figure(2); clf;
plot(b,h)

%% now use the rest of them
num_evos = length(unique(xv(four_inds)));
%fdata = squeeze(mean(reshape(four_data,size(four_data,1),length(four_inds)/num_evos,num_evos),2));
fdata = squeeze(mean(reshape(four_data,size(four_data,1),num_evos, length(four_inds)/num_evos),3));
% for j = 1:size(data,1)
%    fdata(j,:) = mean(reshape(data(j,four_inds),num_evos,length(four_inds)/num_evos),2); 
% end
evotime = mean(reshape(xv(four_inds),num_evos,length(four_inds)/num_evos),2);

L = length(evotime);
T_samp = abs(diff(evotime(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(L);
ffreqs = (dbz_avg/(2*pi))+(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
%prior = exp(-(ffreqs-omegas).^2/(2*new_ddbz^2));prior = prior/sum(prior);


for j = 1:size(fdata,1)
    if 1
    ft = fft(fdata(j,:)-mean(fdata(j,:)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    %figure(4); clf; plot(ffreqs,ft); pause
    [~,mi]=max(ft);
    %omegas_new(j) = omegas(j)+ffreqs(mi);
    newscale(j) = ffreqs(mi);
    else
        Bfreqs = omegas(j)+new_ddbz*linspace(-2,2,1000);
        prior = exp(-(Bfreqs-omegas(j)).^2/(2*new_ddbz^2));prior = prior/sum(prior);
        P_MLE = prior; 
        %figure(4); clf; plot(Bfreqs,P_MLE); pause
        for k = 1:size(fdata,2)
        P_MLE = P_MLE.*(1+fdata(j,k)*cos(Bfreqs*evotime(k)))/(sum(P_MLE));
        %figure(4); clf; plot(Bfreqs,P_MLE); pause
        end
        [~,mi]=max(P_MLE);
        newscale(j) = Bfreqs(mi);
        newvar(j) = std(P_MLE);
    end
end
newscale = -.5+(2*pi*newscale)/dbz_avg;

%%
%newscale = newscale./omegas; %newscale = newscale-mean(newscale);
%newscale = newscale-mean(newscale); newscale = newscale/range(newscale);
%alpha = linspace(-.5,.5,200);
clear tsV ts dataI;
clear dBx
tfunc = @(t,a,s) t*(1/((a+s+.5))); alpha = .5*linspace(-1,1,500);
for i=1:length(alpha);
    for j=1:Nrep;
        %ts=t/(1+alpha(i)*newscale(j));
        %ts=t*(1+alphabest*scale(j))*(1+alpha(i)*newscale(j));
        %ts=t*(1+alpha(i)*dbz_avg/newscale(j));
        ts = tfunc(t,alpha(i),newscale(j));
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=mean(dataI);
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);

cc=[];
rrng = 2:size(talpha,1);
beta0 = [.3, .3, 14, .001, 80 0 1.1];
for j = 1:length(alpha)
    mask = ~isnan(dBx(rrng,j));
    %beta0 = [.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(mask,j))',dBx(mask,j)',beta0,fitfn, [1 1 1 1 1 0]);
    beta0 = pars;
    cc(j) = abs(pars(5));
    %pause
end

cc(cc>1000)=0;
figure(3); clf; plot(alpha,cc);


[T2,mi]=max(cc);
alphabest = alpha(mi);



%%
num_evos = length(unique(xv(shot_inds)));
%fdata = squeeze(mean(reshape(four_data,size(four_data,1),length(four_inds)/num_evos,num_evos),2));
fdata2 = squeeze(mean(reshape(shot_data,size(shot_data,1),num_evos, length(shot_inds)/num_evos),3));
% for j = 1:size(data,1)
%    fdata(j,:) = mean(reshape(data(j,four_inds),num_evos,length(four_inds)/num_evos),2); 
% end

evotime = mean(reshape(xv(shot_inds),num_evos,length(shot_inds)/num_evos),2);

L = length(evotime);
T_samp = abs(diff(evotime(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(L);
ffreqs2 = (F_samp/2)*linspace(0,1,NFFT/2+1);
%prior = exp(-(ffreqs-omegas).^2/(2*new_ddbz^2));prior = prior/sum(prior);


for j = 1:size(fdata2,1)
    if 0
    ft = fft(fdata2(j,:)-mean(fdata2(j,:)),NFFT)/L;
    ft = 2*abs(ft(1:NFFT/2+1));%.*(prior);
    figure(4); clf; plot(ffreqs2,ft); pause
    [~,mi]=max(ft);
    %omegas_new(j) = omegas(j)+ffreqs(mi);
    newscale2(j) = ffreqs2(mi);%+ (newscale(j)+.5)*dbz_avg/(2*pi);
    else
        Bfreqs = omegas(j)+new_ddbz*linspace(-2,2,1000);
        prior = exp(-(Bfreqs-omegas(j)).^2/(2*new_ddbz^2));prior = prior/sum(prior);
        P_MLE = prior; 
        %figure(4); clf; plot(Bfreqs,P_MLE); pause
        for k = 1:size(fdata2,2)
        P_MLE = P_MLE.*(1+fdata(j,k)*cos(Bfreqs*evotime(k)))/(sum(P_MLE));
        %figure(4); clf; plot(Bfreqs,P_MLE); pause
        end
        [~,mi]=max(P_MLE);
        newscale2(j) = Bfreqs(mi);
        newvar(j) = std(P_MLE);
    end
end
newscale3 = -.5+(newscale2)/dbz_avg;

%%
tfunc = @(t,a,s) t*(1/((a+s+.5))); alpha2 = .5*linspace(-.1,4,500);
for i=1:length(alpha2);
    for j=1:Nrep;
        %ts=t/(1+alpha(i)*newscale(j));
        %ts=t*(1+alphabest*scale(j))*(1+alpha(i)*newscale(j));
        %ts=t*(1+alpha(i)*dbz_avg/newscale(j));
        ts = tfunc(tfunc(t,alphabest,newscale(j)),alpha2(i),newscale3(j));
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=mean(dataI);
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha2,ti);

ccc=[];
rrng = 2:size(talpha,1);
beta0 = [.3, .3, 14, .001, 80 0 1.1];
for j = 1:length(alpha2)
    mask = ~isnan(dBx(rrng,j));
    %beta0 = [.3, .3, 14, .001, 80 0 1.1];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(mask,j))',dBx(mask,j)',beta0,fitfn, [1 1 1 1 1 0]);
    beta0 = pars;
    ccc(j) = abs(pars(5));
    newQ(j) = abs(pars(5)/pars(3));
    %pause
end

ccc(ccc>1000)=0;
newQ(newQ>50)=0;
figure(4); clf;
plotyy(alpha2,ccc,alpha2,newQ);
figure(3); clf; plot(alpha2,ccc);


[T2new,mi]=max(cc);
alpha2best = alpha2(mi);
%%
%%
%%
w0 = .082*2*pi;
t = 20:20:20*50;%20:20:50*20;
s = cos(w0*t)+0.2*randn(1,length(t));
t2 = randperm(1000,250);
s2 = cos(w0*t2)+.002*randn(1,length(t2));

L = length(s);
T_samp = abs(diff(t(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(L);
ffreqs = (w0/(2*pi))+(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
%ffreqs = (F_samp/2)*linspace(0,1,NFFT/2+1);
%ffreqs = (F_samp/2)*linspace(-1,1,length(s));

%ft = abs(fftshift(fft(s)));
ft = fft(s); 
ft = 2*abs(ft(1:NFFT/2+1)); ft(end)=0;

figure(1); clf;
plot(t,s);
figure(2); clf;
plot(ffreqs,ft); hold on;
YL = get(gca,'YLim');
[~,mi]=max(ft);
nyquist_limit_freq = ffreqs(mi);
plot(nyquist_limit_freq*[1 1],YL,'r');
fprintf('quantization error = %f MHz\n',1e3*abs(nyquist_limit_freq-w0/(2*pi)))

fine_freqs = linspace(0,ffreqs(end),20*length(ffreqs));


%%
A = zeros(length(ffreqs),length(fine_freqs));
tic
for j = 2:length(ffreqs)
    for k = 2:length(fine_freqs)
        x = linspace(0,2*pi/fine_freqs(2),10^3);
        y = cos(ffreqs(j)*x).*cos(fine_freqs(k)*x);
        A(j,k) = trapz(x,y);
    end  
end
A = A*1/pi/length(x);
A(1,1) = 1;
toc

%%
opts = spgSetParms('verbosity',0);
x = spg_bp(A, ft', opts);
[~,mi_new]=max(abs(x));
new_freq = fine_freqs(mi_new);
figure(3); clf; hold on;
plot(fine_freqs,abs(x)/max(abs(x)),'b','DisplayName','Compressive Sensing');
plot(ffreqs,ft/max(ft),'r','DisplayName','Nyquist limit sensing');
YL = get(gca,'YLim');
plot([1 1]*w0/(2*pi),YL,'k','DisplayName','Real \Omega');
%legend show
