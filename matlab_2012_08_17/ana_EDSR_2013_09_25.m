%%
files = get_files('sm*.mat');

for j = 1:length(files)
d=ana_avg(files{j}); close([1 11 401]); 
figure(200+j); clf
imagesc(d.xv{1},1e-9*linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints),squeeze(nanmean(d.data{1}))); 
ylabel('Frequency (GHz)'); xlabel('epsilon (mV)');
d2=load(d.filename,'configvals');
[T,P,R]=cart2sph(d2.configvals(end-2),d2.configvals(end-1),d2.configvals(end));
%title(sprintf('%.2f T phi = %.0f^o, theta = %.0f^o',[R,90-(P*180/pi),T*180/pi]))
pp=plsinfo('params',d.scan.data.pulsegroups.name,[],d.scantime);
title(sprintf('%.2f T phi = %.0f^o, theta = %.0f^o \n wait time = %.0fns',[R,90-(P*180/pi),T*180/pi,pp(2)]))
end

%%
files = get_files('sm*.mat');
c='rgbcmyk'; c = [c c c c c];
d= ana_avg(files); 
close([1 11 401]); close all
clear h; clear h2;
for j = 1:length(files)
d2=load(d(j).filename,'configvals');
[T,P,R]=cart2sph(d2.configvals(end-2),d2.configvals(end-1),d2.configvals(end));
d(j).phi = 90-(P*180/pi);
h(j) = max(abs(squeeze(nanmean(d(j).data{1}))));
h2(j) = max(abs(squeeze(nanmean(d(j).data{1}(:,:,56:end)))));
end
[~, II] = sort([d.phi]);
figure(400); clf; hold on;
for j = II
    plot(squeeze(nanmean(d(j).data{1})),c(j),'DisplayName',sprintf('phi = %.0f degrees',d(j).phi));
end
legend show;
figure(401);
plot([d.phi],h,'x','MarkerSize',10,'LineWidth',3);
xlabel('\phi');
ylabel('height of ST+ peak');

figure(402);
plot([d.phi],h2,'x','MarkerSize',10,'LineWidth',3);
xlabel('\phi');
ylabel('height of second peak');

%%
files = get_files('sm*.mat');

for j = 1:length(files)
d=ana_avg(files{j}); close([1 11 401]); 
figure(300+j); clf
d2=load(d.filename,'configvals');
cntrfrq = d2.configvals(28);
imagesc(1e-9*cntrfrq-(4:4:128*4)*1e-3,d.tv,squeeze(nanmean(d.data{1}))); 
set(gca,'YDir','normal')
set(gca,'XDir','reverse')
xlabel('Frequency (GHz)'); ylabel('epsilon (mV)');
[T,P,R]=cart2sph(d2.configvals(end-2),d2.configvals(end-1),d2.configvals(end));
title(sprintf('%.2f T phi = %.0f^o, theta = %.0f^o',[R,90-(P*180/pi),T*180/pi]))
%pp=plsinfo('params',d.scan.data.pulsegroups.name,[],d.scantime);
%title(sprintf('%.2f T phi = %.0f^o, theta = %.0f^o \n wait time = %.0fns',[R,90-(P*180/pi),T*180/pi,pp(2)]))
end

%%
% Use for extracting line cuts from data sets with varying ARP time length
files = get_files('sm*.mat');
lines=0;
cc = 'rgbcmyk'; cc=[cc cc cc cc];
figure(400); clf;
for j = 1:length(files)
opt='';
%opt='noscale';
d=ana_avg(files{j},opt);
try
    close(1);
    close(11);
    %close(401);
end
figure(200+j); clf
imagesc(d.xv{1},1e-9*linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints),squeeze(nanmean(d.data{1}))); 
colorbar;
%set(gca,'YDir','normal');
datatemp=squeeze(nanmean(d.data{1}));
slice=4.956e9;
increment=(d.scan.loops(1).rng(2)-d.scan.loops(1).rng(1))/d.scan.loops(1).npoints;
index=round((slice-d.scan.loops(1).rng(1))/increment);
lines(j,1:64)=datatemp(index,1:64);
ylabel('Frequency (GHz)'); xlabel('ARP time (us)');
fnames{j}=d.filename;
filename=d.filename;
filename(strfind(filename,'_'))='';
xdata=d.xv{1};
ydata=lines(j,:);
d2=load(d.filename,'configvals');
Pow=d2.configvals(end-3);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
eps1=s(5,1);
eps2=s(4,1);
%title(sprintf('%.2f T phi = %.0f^o, theta = %.0f^o',[R,90-(P*180/pi),T*180/pi]))
pp=plsinfo('params',d.scan.data.pulsegroups.name,[],d.scantime);
title(sprintf('%s Power = %.2f dBm \n eps1 = %.2f eps2 = %.2f %s',filename,Pow,eps1,eps2,opt))
figure(400);
hold on;
plot(d.xv{1},lines(j,1:64),cc(j),'DisplayName',sprintf('%.2f dBm eps1= %.2f mV eps2 = %.2f mV',Pow,eps1,eps2));
xlabel('ARP time (us)');
ylabel('Transition prob')
title(sprintf('Transition probability at %.2f GHz %s',slice/1e9,opt));
%leg{j}=
%legend(leg);
end
legend show

%% Just read the files and graph the data

files = get_files('sm*.mat')
for j=1;length(files)
%opt='';
opt='noscale';
d=ana_avg(files{j},opt); 
try
    close(1);
    close(11);
    close(401);
end
figs(j)=200+j;
figure(200+j); clf
imagesc(d.xv{1},1e-9*linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints),squeeze(nanmean(d.data{1}))); 
colorbar;
%set(gca,'YDir','normal');
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
eps1=s(5,1);
eps2=s(4,1);
ylabel('Frequency (Hz)'); xlabel('ARP time (us)');
fnames{j}=d.filename;
filename=d.filename;
filename(strfind(filename,'_'))='';
title(sprintf('%s \n eps1 = %.1f mV eps2 =  %.1f mV %s',filename,eps1,eps2,opt));
ppt=guidata(pptplot);
%set(ppt.e_file,'String',filename{1});
set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
%set(ppt.e_title,'String',s(1).scan.data.prettyname);
%set(ppt.e_body,'String','');
%set(ppt.exported,'Value',0);
end


%% spyview some EDSR scans eps vs freq

d=ana_avg;
s_1=spyview(squeeze(nanmean(d.data{1})),d.xv{1},1e-9*linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints));

%% now do the rest
fts = s_1.filters;
close all;
files = get_files;
clear s;
for j = 1:length(files)
    d = ana_avg(files{j});
    s(j)=spyview(squeeze(nanmean(d.data{1})),d.xv{1},1e-9*linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints));
    s(j).filters = fts;
end
