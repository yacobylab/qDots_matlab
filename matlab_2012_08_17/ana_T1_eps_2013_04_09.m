%% get the files
files = uigetfile('sm*.mat','MultiSelect','on');
if ~iscell(files) 
    files = {files};
end
files = sort(files);
% {sm_T1_eps_R_3582.mat,
% sm_T1_eps_R_3583.mat,
% sm_T1_eps_R_3584.mat,
% sm_T1_eps_R_3585.mat,
% sm_T1_eps_R_3586.mat,
% sm_T1_eps_R_3587.mat,
% sm_T1_eps_R_3588.mat,
% sm_T1_eps_R_3589.mat,
% sm_T1_eps_R_3590.mat};
%% analyze them
fitfn = @(p,x) p(1)+p(2)*(1-exp(-1*x/p(3)));
results = [];
for j = 1:length(files)
s=ana_avg(files{j}); close all;
scantime = getscantime(s.scan,s.data);
data = squeeze(nanmean(s.data{1}));
gs = data(2,:);
pinfo = plsinfo('xval',s.scan.data.pulsegroups(2).name,[],scantime);
ts = pinfo(2,:);
eps = pinfo(1,1);
figure(1); clf;
plot(ts,gs); hold on;
%pars = fitwrap('plinit plfit',ts,gs,[.01, range(gs), ts(end/2)],fitfn);
pars = polyfit(ts,gs,1);
plot(ts,polyval(pars,ts),'r');
results(end+1).T1 = 1/pars(1);

% if pars(3)>1000 || pars(3) < 0
%    results(end).T1 = 0; 
% end
results(end).eps = eps;
results(end).ts = ts;
results(end).data = gs;
results(end).file = files{j};
end

%% Set up J(eps) and the related

%sm_Ramsey_R_3581.mat
[a b c]=ana_echo('',struct('opts','ramsey','grng', [2.1 4]));
djde = @(eps) b.freqfunc(eps)/b.decp(2);
dbz = 1/(1e-9*2*pi/b.omega_dbz); %Hz;
%% Plot Nice stuff and extract noise

js = 1e6*b.freqfunc([results.eps]); %Hz
djvals = 1e9*djde([results.eps]); %Hz/V
omegas = 1e6*sqrt(js.^2+dbz^2); %Hz
seps = 1./(1e-6*[results.T1].*djvals.^2.*dbz./sqrt(omegas));
figure(1); clf; 
logj = log(js); 
logd=log(seps);
plot(logj,logd,'bx'); hold on;
noise = polyfit(logj,logd,1);
plot(logj,polyval(noise,logj),'r');
xlabel('Log(J)'); ylabel('Log S_\epsilon');
OneMHzNoise = sqrt(exp(polyval(noise,log(1e6)))); 
title(sprintf('Noise at 1MHz = %.2f nV/sqrt(Hz)\n noise ~ 1/f^{%.1f}',1e9*OneMHzNoise,-noise(1)));
figure(2); clf;
plot(js,seps,'bx','MarkerSize',7);
hold on;
xlabel('J(Hz)'); ylabel('S_\epsilon (V^2/Hz)');
plot(js,exp(polyval(noise,logj)),'r');
title(sprintf('Noise at 1MHz = %.2f nV/sqrt(Hz)\n noise ~ 1/f^{%.1f}',1e9*OneMHzNoise,-noise(1)));


