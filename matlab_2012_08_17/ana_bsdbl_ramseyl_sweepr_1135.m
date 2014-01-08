%
load(uigetfile('sm_bsdbl*.mat'));
data=data{1};
data=squeeze(nanmean(data,1));

T=(1:101);
mask=12:52;
fitfn=@(p,x) p(1).*exp(-x.^2./p(2)).*cos(p(3).*x+p(4))+p(5);
beta0=[-1e-3 1e3 1 0 -8.5e-3];
params=[];
paramtmp=[];
T=T(mask);
for i=1:size(data,1)
   
datafit=data(i,mask);
guess=fioscill(T,datafit,2);
beta0(3)=guess(4);
paramtmp(i,:) = fitwrap('plinit plfit',T,datafit, beta0, fitfn, [1 1 1 1 1]);
end

    paramtmp=paramtmp./(2*pi).*1000;
figure(200)
x=linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),size(paramtmp,1));
plot(x*1e3,paramtmp(:,3))
xlabel('Epsilon (mv)')
ylabel('Frequency (MHz)')

%%
figure(200)
clf
offset=18e-4;
i=6
datatmp=data(i,mask)
  plot(T,datatmp,'b.-')
    hold on
i=46
datatmp=data(i,mask)+offset
  plot(T,datatmp,'m.-')


i=

%% waterfall 
figure(400)
clf
offset=0;
colors=[{'g'},{'m'},{'c'},{'r'},{'b'},{'y'}]
for i=1:2:size(data,1)
    cind=mod((i+1)/2,6)+1
    offset=offset+12e-4
    datatmp=data(i,mask)+offset
    plot(T,datatmp,colors{cind})
    hold on
end

xlabel('Time (ns)')
ylabel('Epsilon (mV)')
ylim=([scan.loops(1).rng(1),scan.loops(1).rng(2)]);
    
%%
load(uigetfile('sm_bsdbl*.mat'));
data=data{1};
data=squeeze(nanmean(data,1));
T=(1:101);
mask=15:50;
fitfn=@(p,x) p(1).*exp(-x.^2./p(2)).*cos(p(3).*x+p(4))+p(5);
beta0=[-1.5e-3 1e3 1 0 -7.5e-3];
params=[];
paramtmp=[];
T=T(mask);
for i=1:size(data,1)
   
datafit=data(i,mask);
guess=fioscill(T,datafit,2);
beta0(3)=guess(4);
paramtmp(i,:) = fitwrap('plinit plfit',T,datafit, beta0, fitfn, [1 1 1 1 1]);
end

    paramtmp=paramtmp./(2*pi).*1000;
figure(200)
x=linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),size(paramtmp,1));
plot(x*1e3,paramtmp(:,3))
xlabel('Epsilon (mv)')
ylabel('Frequency (MHz)')
