

%% Graph data 
% use for jumping and holding
files = get_files('sm*.mat')
for j=1;length(files)
%d=ana_avg(files{j},'noscale'); 
d=ana_avg(files{j}); 
try
    close(1);
    close(11);
    %close(401);
end
data=squeeze(nanmean(d.data{1}));
% data2=squeeze(nanmean(d2.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
% s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);
% eps2s=s(4,1);
% eps2f=s2(4,1);
% eps1=s2(5,1);
eps2s=s(1,1);
eps2f=s2(1,1);
eps1=s2(2,1);
data=data(:,6:end);
yvals=linspace(eps2s,eps2f,ngroups);
xvals=d.xv{1};
figure(200+j); clf
% subplot(1,2,1);
imagesc(xvals,yvals,data); 
colorbar;
xlabel('time (ns)');
ylabel(sprintf('eps (mV',eps1));
% eps=(eps1);
% slice=hslice(200+j,eps);
% %figure(200+j);
% %subplot(1,2,2);
% plot(slice(1,:),slice(2,:))

title(sprintf('Linecut at eps = %.2f',eps));
xlabel('time (ns)');
B=data-mean(mean(data));
size(B)
B=fft(B');
figure(270+j);
imagesc(real(B'));
caxis([-.5 .5]);
figure(270+j+1);
imagesc(imag(B'));
caxis([-.5,.5]);
% figure(250+j); clf
% imagesc(xvals,yvals,data2); 
% colorbar;
% xlabel('time (us)');
% ylabel(sprintf('eps2 (mV) (eps1 = %.2f)',eps1));
end


%% Graph, fft, and reject early times
% use for jumping and holding
files = get_files('sm*.mat')
for j=1;length(files)
%d=ana_avg(files{j},'noscale'); 
d=ana_avg(files{j}); 
try
    close(1);
    close(11);
    close(401);
end
data=squeeze(nanmean(d.data{1}));
figure(1);
imagesc(data);
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);
eps2s=s(1,1);
eps2f=s2(1,1);
%take only the data where the pulse has turned on. For each row, loop
%through the columns and find where it turns on. Take only the data
%afterward and zero pad the rest. 

start=1;
dataf=data.*0;
threshold=0;
for i=1:size(dataf,1)
    if ~isempty(find(abs(data(i,1:end))>threshold,1))
       ind(i)=find(abs(data(i,1:end))>threshold,1);
    else
      ind(i)=1;
    end
    
    dataf(i,1:end-ind(i)+1)=data(i,ind(i):end);
end

%put a window on the data
win=exp(-([1:size(dataf,2)]/40).^2);
win=repmat(win,size(dataf,1),1);
dataf=diff(data')';
% dataf=dataf.*win;
% dataf=dataf(start:end,:);

% for j=1:size(data,1);
%     p=polyfit(1:1:64,data(j,:),4);
%     datap(j,:)=polyval(p,1:1:64);
% end
% dataf=data-datap;


yvals=linspace(eps2s,eps2f,ngroups);
xvals=d.xv{1}*1;
figure(201); clf
imagesc(xvals,yvals,data); 
colorbar;
xlabel('time (ns)');
ylabel('eps (mV)');

figure(221);
imagesc(xvals,yvals,dataf); 
colorbar;
xlabel('time (ns)');
ylabel('eps (mV)');
A=data;
B=dataf;
size(B);
A=fft(A');
B=fft(B');
% figure(261);
% imagesc(real(B'));
% title('Real');
% colorbar;
% caxis([-.5 .5]);
figure(262);
imagesc(abs(A'));
title('raw Magnitude');
colorbar;
caxis([0,1]);
figure(263);
imagesc(abs(B'));
title('filtered Magnitude');
caxis([0 1]);
colorbar;
montagefigures([201 221 262 263], 2,2);

end
%% Graph data for arps
files = get_files('sm*.mat')
for j=1;length(files)
%d2=ana_avg(files{j},'noscale'); 
d=ana_avg(files{j}); 
try
    close(1);
    close(11);
    %close(401);
end
xvp=6;
xlab='time (us)';
yvp=4;
ylab='eps (mV)';
data=squeeze(nanmean(d.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);

yvals=linspace(s(yvp,1),s2(yvp,1),ngroups);
xvals=d.xv{1};
figure(j); clf
imagesc(xvals,yvals,data); 
colorbar;
xlabel(xlab);
ylabel(ylab);

end

%% 
files = get_files('sm*.mat')
for j=1;length(files)
opt='';
%opt='noscale';
d=ana_avg(files{j},opt); 
try
    close(1);
    close(11);
    close(401);
end
data=squeeze((d.data{1}));
%yvals=linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
yvals=linspace(1,size(squeeze(d.data{1}),1));
xvals=d.xv{1};
figure(300+j); clf
imagesc(xvals,yvals,data); 
colorbar;
a=gca;
d=caxis(a);
d(1)=0;
caxis(d);
set(gca,'YDir','normal');
xlabel('epsilon (mV)');
ylabel('count');
%vslice(300+j,.811,0.06);
end

%% 
files = get_files('sm*.mat')
for j=1;length(files)
opt='';
%opt='noscale';
d=ana_avg(files{j},opt); 
try
    close(1);
    close(11);
    close(401);
end
data=squeeze(nanmean(d.data{1}));
yvals=3*1e9./d.scan.data.freqs;
xvals=d.xv{1};
figure(300+j); clf
imagesc(xvals,yvals,data); 
colorbar;
xlabel('epsilon (mV)');
ylabel('evolution time (ns)');
vslice(300+j,.84,0.05);
end

%% Graph data for freq vs epsilon sweeps.
opt='';
%opt='noscale';
d=ana_avg(opt); 
try
    close(1);
    close(11);
    close(401);
end
data=squeeze(nanmean(d.data{1}));
yvals=linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
xvals=d.xv{1};
figure(401); clf
imagesc(xvals,yvals,diff(data')');
figure(501); clf
imagesc(xvals,yvals,diff(data));

ngroups=size(d.scan.data.pulsegroups,2);
gd=plsinfo('gd',d.scan.data.pulsegroups(ceil(ngroups/2)).name,[],d.scantime);
output=atplschk2(gd.name,size(gd.varpar,1)/2);
outchans=output(:,1);
outmark=output(:,2);
close(55);


figure(301); clf
hold on;
subplot(3,1,[1 2])
imagesc(xvals,yvals,data); 
colorbar;
xlabel('epsilon (mV)');
ylabel('Freq');
subplot(3,1,3);

colors={[1 0 0],[0.5 0 0], [0 1 0], [0 0.5 0]}; 
colorsm={[1 0 0],[.25 0 0],[1 0 0],[1 1 0]};
for l=1:length(outchans)
p=plot(outchans{l});
hold on;
pm=plot(outmark{l});
set(p,'Color',colors{l});
set(pm,'Color',colors{l});
end

%% Graph data for freq vs epsilon sweeps.
%subtract a background file.
opt='';
%opt='noscale';
d=ana_avg(opt); 
try
    close(1);
    close(11);
    close(401);
end

b=ana_avg(opt);
try
    close(1);
    close(11);
    close(401);
end
b=squeeze(nanmean(b.data{1}));
bavg=squeeze(nanmean(b));
b=repmat(bavg,size(b,1),1)

data=squeeze(nanmean(d.data{1}))-b;
yvals=linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
xvals=d.xv{1};
figure(401); clf
imagesc(xvals,yvals,data);

ngroups=size(d.scan.data.pulsegroups,2);
gd=plsinfo('gd',d.scan.data.pulsegroups(ceil(ngroups/2)).name,[],d.scantime);
output=atplschk2(gd.name,size(gd.varpar,1)/2);
outchans=output(:,1);
outmark=output(:,2);
close(55);


figure(301); clf
hold on;
subplot(3,1,[1 2])
imagesc(xvals,yvals,data); 
colorbar;
xlabel('epsilon (mV)');
ylabel('Freq');
subplot(3,1,3);

colors={[1 0 0],[0.5 0 0], [0 1 0], [0 0.5 0]}; 
for l=1:length(outchans)
p=plot(outchans{l});
hold on;
pm=plot(outmark{l});
set(p,'Color',colors{l});
set(pm,'Color',colors{l});
end


%% Graph the pulse group for a specific dataset.
files = get_files('sm*.mat')
for j=1;length(files)
d=ana_avg(files{j}); 
try
    close(1);
    close(11);
    close(401);
end
data=squeeze(nanmean(d.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
gd=plsinfo('gd',d.scan.data.pulsegroups(ceil(ngroups/2)).name,[],d.scantime);
figure(55); clf;
atplschk2(gd.name,ceil(size(gd.varpar,1)/2));
end

%% General display data for multiple pulsegroups
files = get_files('sm*.mat')
for j=1;length(files)
%d2=ana_avg(files{j},'noscale'); 
d=ana_avg(files{j}); 
try
    close(1);
    close(11);
    %close(401);
end
xvp=2;
xlab='time (us)';
yvp=1;
ylab='eps (mV)';
data=squeeze(nanmean(d.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);

yvals=linspace(s(yvp,1),s2(yvp,1),ngroups);
xvals=d.xv{1};
figure(j); clf
imagesc(xvals,yvals,data); 
colorbar;
xlabel(xlab);
ylabel(ylab);

end



%% Analyze TpT1 data

files = get_files('sm*.mat')
j=1;
d=ana_avg(files{j}); 
try
    %close(1);
    close(11);
    %close(401);
end
xvp=2;
xlab='time (us)';
yvp=1;
ylab='eps (mV)';
data=squeeze(nanmean(d.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);
yvals=linspace(s(yvp,1),s2(yvp,1),ngroups);
xvals=d.xv{1};
figure(1);clf;
imagesc(xvals,yvals,data);

datasub=data(yvals<1.0 & yvals>0,:);
ysub=yvals(yvals<1.0 & yvals>0);

clear tau;
clear amp;
clear offset;
clear fitdata;
for j=1:size(ysub,2)
    beta0=[.2 .8 4];
    slice=datasub(j,2:end);
    xslice=xvals(2:end);
    fitfunc = @(b,x) b(1)+b(2)*exp(-x/b(3));
    beta=nlinfit(xslice,slice,fitfunc,beta0);
        figure(2); clf;
    plot(xslice,slice);
    hold on;
    fit=fitfunc(beta,xslice);
    plot(xslice,fit);
    fitdata(j,:)=fit;
    figure(3);clf;
    imagesc(fitdata);
    tau(j)=beta(3);
    amp(j)=beta(2);
    offset(j)=beta(1);
end


%% Analyze adiabatic passage data to extract the size of the STp splitting.

%887 mT conversion
fx1='sm_STP_R_8393';
fz1='sm_STP_R_8408';%This is actually at 700 mT
c1=.07; %mV/GHz

%400 mT conversion
fx2='sm_STP_R_8397';
fz2='sm_STP_R_8406';
c2=.14; %mV/GHz

%100 mT conversion
fx3='sm_STP_R_8400';
fz3='sm_STP_R_8403';
c3=15;

files={fx1,fx2,fx3,fz1,fz2,fz3};
fields=[887 400 100 700 400 100];
c=[c1 c2 c3 c2 c2 c3]*1e-9;

xepsmin=[1.4 1.6 2];
xepsmax=[2.2 3 3];

zepsmin=[1.5 1.5 1.5];
zepsmax=[2.5 2.5 2.5];

epsmin=[xepsmin zepsmin];
epsmax=[xepsmax zepsmax];

xvp=6;
yvp=4;

for i=1:6
d=ana_avg(files{i});close(1);close(11);close(401);
data=squeeze(nanmean(d.data{1}));
ngroups=size(d.scan.data.pulsegroups,2);
s=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
s2=plsinfo('xval',d.scan.data.pulsegroups(ngroups).name,[],d.scantime);
yvals=linspace(s(yvp,1),s2(yvp,1),ngroups);
xvals=d.xv{1};


datasub=data(yvals<epsmax(i) & yvals>epsmin(i),:);
ysub=yvals(yvals<epsmax(i) & yvals>epsmin(i));

for j=1:size(ysub,2)
    beta0=[.5 .5 .1];
    slice=datasub(j,:);
    fitfunc = @(b,x) b(1)-b(2)*exp(-x/b(3));
    beta=nlinfit(xvals,slice,fitfunc,beta0);
    timeconsts(j)=beta(3)*1e-6;
    rabifreqs(j)=sqrt(ysub(j)/(pi^2*timeconsts(j)*c(i)));
end

splitting(i)=mean(rabifreqs);
splitstdev(i)=std(rabifreqs)/sqrt(size(rabifreqs,2));

end

xsplit=splitting(1:3);
xfields=fields(1:3);
xerror=splitstdev(4:6);
zsplit=splitting(4:6);
zfields=fields(4:6);
zerror=splitstdev(4:6);

figure(1);clf;
errorbar(xfields,xsplit,xerror,'r','DisplayName','X field');
hold on;
errorbar(zfields,zsplit,zerror,'b','DisplayName','Z field');
xlabel('Field (mT)');
ylabel('Splitting Hz');
legend show;






%% Compute Rabi rate
d_omegad_e=2*pi*5e9/(.6);
d_eps=2;
d_omega=d_omegad_e*d_eps;
d_t=100e-9;
%Rabi frequency
fR=sqrt(d_omega*2/(d_t*pi))/(2*pi);
% Rabi Period
TR=1/fR;
% Energy in eV
E=fR*2*pi*6.6e-34/(2*pi)/(1.6e-19)

%% Compute splitting based on fitting the data at the stp point
gamma=2*pi*5.6e9;
T2=(1/.06887)*1e-9;
Bnuc=1/(T2*gamma)
h=6.6e-34;
hbar=h/(2*pi);
E=gamma*Bnuc*hbar/(1.6e-19)

%% 
Rseries=247;
Rshunt=61;
Rload=50;
Requiv=(1/Rshunt+1/((1/Rload+1/Rshunt)^(-1)+Rseries))^(-1)




