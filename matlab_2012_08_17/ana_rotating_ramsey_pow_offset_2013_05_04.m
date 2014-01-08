
%%
files = uigetfile('sm*amp*.mat','MultiSelect','on');
if ~iscell(files)
    files = {files};
end
files = sort(files);

%%
results = [];

for i = 1:length(files)
d=ana_avg(files{i},'noppt');
pp=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
xv = pp(2,:);
results(i).eps_off = pp(1,1);
data = squeeze(nanmean(d.data{1}));
results(i).data = data;
powvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);

%d= ana_avg(file,'noppt');
%data = squeeze(nanmean(d.data{1}));
%fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
freq = []; ampl=[]; t2s = [];
%pp=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
%xv = pp(2,:);
%fvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),size(data,1));

cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
fifn.fn = @fioscill;
fifn.args = {2};
mask = 12:size(data,2);
for j = 1:size(data,1)
     fifn.args={2};
     fp=fitwrap('fine woff',xv(mask),data(j,mask),fifn, cosfn5, [1 0 1 1 0 0]);
     fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
     pars=fitwrap('plinit plfit fine woff',xv(mask),data(j,mask),fp, cosfn2, [1 1 1 1 0 1]);
     freq(j) = abs(1e9*pars(4)/(2*pi)); %abs(1/pars(3));
     ampl(j) = sqrt(pars(2)^2+pars(3)^2);%abs(pars(2));
     t2s(j) = 1e-9/pars(6);
end
results(i).freq = freq;
results(i).ampl = ampl;
results(i).t2s = t2s;
results(i).filename = files{i};
results(i).pars = pars;
end
%%
close all;
allfreqs = (vertcat(results.freq));
allt2s = abs(vertcat(results.t2s));
allfreqs(allfreqs>70e6) = NaN;
allt2s(allt2s>700e-9) = NaN;
Qs = allfreqs.*allt2s;
Qs(Qs>15) = NaN;

figure(1); clf;
imagesc([results.eps_off],powvals,1e-6*allfreqs(end:-1:1,end:-1:1));
caxis([5 50]);
xlabel('\epsilon_{offset}');
ylabel('RF power (dBm from generator)');
title('Rotating Ramsey Frequency (MHz)');
colorbar

figure(2); clf;
imagesc([results.eps_off],powvals,1e9*allt2s(end:-1:1,end:-1:1));
caxis([5 550]);
xlabel('\epsilon_{offset}');
ylabel('RF power (dBm from generator)');
title('T2^* (ns)');
colorbar

figure(3); clf;
imagesc([results.eps_off],powvals,Qs(end:-1:1,end:-1:1));
caxis([2 14]);
xlabel('\epsilon_{offset}');
ylabel('RF power (dBm from generator)');
title('Q');
colorbar

%%
cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
data = results(1).data; mask  = 6:size(data,2);
for j = 1:size(data,1)
     fifn.args={2};
     fp=fitwrap('fine woff',xv(mask),data(j,mask),fifn, cosfn5, [1 0 1 1 0 0]);
     fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
     pars=fitwrap('plinit plfit fine woff',xv(mask),data(j,mask),fp, cosfn2, [1 1 1 1 0 1]);
     freq(j) = abs(1e9*pars(4)/(2*pi)); %abs(1/pars(3));
     ampl(j) = sqrt(pars(2)^2+pars(3)^2);%abs(pars(2));
     t2s(j) = 1e-9/pars(6);
     pause
end