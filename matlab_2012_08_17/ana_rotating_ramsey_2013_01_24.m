%% first frequency sweep
d=ana_avg; close all;
%filename: 'sm_ramsey_rotate_R_3156.mat'
%filename: 'sm_ramsey_rotate_R_3157 .mat'
data = squeeze(nanmean(d.data{1}));
scantime = getscantime(d.scan,d.data);

freq =[];
for j = 2:length(d.scan.data.pulsegroups)
   pp=plsinfo('params',d.scan.data.pulsegroups(j).name,[],scantime);
   freq(end+1) = pp(6);
end
vp=plsinfo('xval',d.scan.data.pulsegroups(j).name,[],scantime);
vp = vp(2,:);


figure(1); clf;
plot(data(1,:));
xlabel('\tau (ns)');
ylabel('P_{Triplet}');
figure(2); clf;

figure(2); clf;
imagesc(vp,freq,data(2:end,:));
xlabel('Burst Time(ns)');
ylabel('Frequency (MHz)');


%% now a power sweep 
d=ana_avg; %close all;
%  filename: 'sm_ramsey_rotate_R_3158.mat'
%  filename: 'sm_ramsey_rotate_R_3159.mat'
%  filename: 'sm_ramsey_rotate_R_3160.mat'
data = squeeze(nanmean(d.data{1}));
scantime = getscantime(d.scan,d.data);
ampl =[];
for j = 2:length(d.scan.data.pulsegroups)
   pp=plsinfo('params',d.scan.data.pulsegroups(j).name,[],scantime);
   ampl(end+1) = pp(5);
end
vp = plsinfo('xval',d.scan.data.pulsegroups(2).name,[],scantime);
vp = vp(2,:);


figure(4); clf;
plot(data(1,:));
xlabel('\tau (ns)');
ylabel('P_{Triplet}');


figure(5); clf;
imagesc(vp,ampl,data(2:end,:));
xlabel('Burst Time(ns)');
ylabel('Amplitude (mV)');


%% more analysis of J(eps)

[a b c]=ana_echo('',struct('opts','ramsey noref','frames',[1:50]));
%close all;
Js = 1e3*b.params(:,4)./(2*pi);
eps = b.freq(1,:);
t2s = abs(1./b.params(:,6))';
mask = find(t2s<.6e3);
%mask=1:length(eps); 
eps = eps(mask);
Js = Js(mask);
t2s = t2s(mask);

figure(1); clf; hold on;
plot(eps, Js,'rx', 'MarkerSize',8, 'LineWidth',2);
xlabel('epsilon (mV)');
ylabel('J (MHz)');
fitfn = @(p,x) p(1)+p(2)*exp(x./p(3));
beta0 = [min(Js), .5*range(Js), .5*range(eps)];
pars = fitwrap('plfit ',eps,Js',beta0,fitfn);
figure(1); 
plot(eps,fitfn(pars,eps),'b');

figure(2); clf; hold on;
plot(eps, t2s,'rx', 'MarkerSize',8, 'LineWidth',2);
xlabel('epsilon (mV)');
ylabel('T_{2}^{*} (ns)');

djfun = @(p,x) (p(2)/p(3))*exp(x./p(3));
djde = djfun(pars,eps);
figure(3); clf; hold on;
plot(1./djde, t2s,'rx', 'MarkerSize',8, 'LineWidth',2);
xlabel('(dj/de)^{-1} (mV/MHz)');
ylabel('T_{2}* (ns)');
noise=polyfit(1./djde,t2s,1); noiseslope = abs(1/noise(1));
plot(1./djde,(noise(1)*(1./djde)+noise(2)),'b')
erms=sqrt(2)*noiseslope/(2*pi);
title(sprintf('T2* = %.0f (dJ/deps)^{-1} + %.0f',noise));
t=text(1/djde(round(end/2)),max(t2s),sprintf('eps_{RMS} = %.2d V',erms));

figure(5); clf; hold on;
plot(eps, 1e-3*t2s.*Js','rx','MarkerSize',8, 'LineWidth',2);
xlabel('epsilon (mV)');
ylabel('Q');

%% generator freq sweep

d = ana_avg;
data = squeeze(nanmean(d.data{1}));
fvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
tvals = plsinfo('xval',d.scan.data.pulsegroups.name,[],d.scantime);
tvals = tvals(2,:);
pinfo = plsinfo('params',d.scan.data.pulsegroups.name,[],d.scantime);
eps = pinfo(4);
figure(102); clf; 
imagesc(tvals,fvals,data); colorbar;
xlabel('evolution time (ns)');
ylabel('Frequency (Hz)');
title(sprintf('epsilon offset = %.4f',eps));

%%
d= ana_avg;
data = squeeze(nanmean(d.data{1}));

%% ypi cal

d=ana_avg; figure(1); clf; imagesc(d.xv{2},d.tv,squeeze(nanmean(d.data{1}(:,2:end,:))));
xlabel('T (ns)');
ylabel('\epsilon (mV)')
