%
d1=load('data_100MHz_samp_11_avg_rf_L.mat');
%d1=load('data_100MHz_samp_avg_norf_L.mat');
%d3=load('data_100MHz_samp_0102_avg_rf_L.mat');
d2=load('data_100MHz_samp_1102_final_avg_L.mat');
d3=load('data_100MHz_samp_1102_3_avg_L.mat');
%d3=load('data_100MHz_samp_0102b_avg_rf_L.mat');
nyquist=50e6; 
figure(1);
sf=20;
clf;
% http://www.mathworks.com/help/signal/ref/dspdata.psd.html for normalizaton
f1=mean(abs(fft(d1.aa,[],2)),1)./sqrt(size(d1.aa,2)*nyquist*2); f1=f1(1:end/2); f1(1)=0;
f2=mean(abs(fft(d2.aa,[],2)),1)./sqrt(size(d1.aa,2)*nyquist*2); f2=f2(1:end/2); f2(1)=0;
f3=mean(abs(fft(d3.aa,[],2)),1)./sqrt(size(d1.aa,2)*nyquist*2); f3=f3(1:end/2); f3(1)=0;
freq=linspace(0,nyquist,length(f1));
loglog(freq,smooth(f1,sf),'DisplayName','11');
hold on;
loglog(freq,smooth(f2,sf),'r-','DisplayName','(11)-(02)');
%loglog(freq,smooth(abs(sqrt(f2.^2-f1.^2)),sf),'r-','DisplayName','Difference');
loglog(freq,smooth(f3,sf),'k-','DisplayName','nognd');
fprintf('Done\n');
legend show;

%%
figure(2);
clf;
nyquist=50e6;
d3=load('data_100MHz_samp_0102_sweep_avg_L.mat');
'loaded'
f3=abs(fft(d3.aa,[],2)); f3=f3(:,1:end/2); f3(:,1)=0;
freq=linspace(0,nyquist,size(f3,2));
'fft'
for i=1:size(f3,1)
    f3s(i,:)=smooth(f3(i,:),80);
end
pcolor(log10(freq),1:size(f3s,1),log(f3s));
shading interp;
figure(3);
plot(mean(d3.aa,2));

%%
figure(1); clf; figure(2); clf; 
smth=80; 
nyquist=50e6; 
colors={'b','r','g','m','c'};
for i =1:4 
    fname=sprintf('sm_100MHz_LR_%d.mat',2110+i);
    data1=load(fname);
    L=length(data1.dd.data{1}); 
    fftl=abs(fft(data1.dd.data{1}))./sqrt(L*2*nyquist); fftl=fftl(1:end/2); fftl(1)=0;
    fftr=abs(fft(data1.dd.data{2}))./sqrt(L*2*nyquist); fftr=fftr(1:end/2); fftr(1)=0;
       note{i}=data1.dd.note;
    freq=linspace(0,nyquist,length(fftl)); 
    figure(1) 
    loglog(freq,smooth(fftl,smth),colors{i}); 
    hold on 
    figure(2) 
    loglog(freq,smooth(fftr,smth),colors{i}); 
    hold on 
end
note{1}='(0,2), sensitive sensors'; note{2}='(0,2), insensitive sensors'; note{3}='(1,1), sensitive sensors'; note{4}='(1,1), insensitive sensors';
figure(1) 
legend(note{1},note{2},note{3},note{4}); 
title('Left Side'); xlabel('Frequency (Hz)'); ylabel('Noise $$(V/\sqrt{Hz})$$','Interpreter','latex','FontSize',13)
figure(2) 
legend(note{1},note{2},note{3},note{4}); 

title('Right Side'); xlabel('Frequency (Hz)'); ylabel('Noise $$(V/\sqrt{Hz})$$','Interpreter','latex','FontSize',13);
    
    
    
%% now compare sides. 
smth=80; 
nyquist=50e6; 
for i =1:4 
    fname=sprintf('sm_100MHz_LR_%d.mat',2110+i);
    data1=load(fname);
    L=length(data1.dd.data{1}); 
    fftl=abs(fft(data1.dd.data{1}))/L; fftl=fftl(1:end/2); fftl(1)=0;
    fftr=abs(fft(data1.dd.data{2}))/L; fftr=fftr(1:end/2); fftr(1)=0;
    freq=linspace(0,nyquist,length(fftl)); 
    figure(i); clf; 
    loglog(freq,1e6.*smooth(fftl,smth),'b'); 
    hold on 
    loglog(freq,1e6.*smooth(fftr,smth),'r'); 
    title(data1.dd.note); xlabel('Frequency (Hz)'); ylabel('Noise (V)');
    legend('left','right')
    
end
%%
figure(5); clf; 
x=1:1e4; 
fftsin=abs(fft(sin(x)));
%fftsin=fftsin(1:end/2); fftsin(1)=0; 
pow=fftsin.^2./1e4;
plot(pow)

%%
figure(5); clf;
r=rand(1,1e4);
fftr=abs(fft(r)); 
fftr(1)=0;
pow=fftr.^2./1e4; 
plot(pow)
%%

[xc, lags] = xcorr(data1.dd.data{1},data1.dd.data{1});
figure(4); clf
plot(xc(1:end/2))

%% 2013_05_18

d1= load('sm_100MHz_LR_4090'); sf=20;
f1=(abs(fft(d1.dd.data{1})))./sqrt(size(d1.dd.data{1},1)*nyquist*2); 
f1=f1(1:end/2); f1(1)=0;
freq=linspace(0,nyquist,length(f1));
figure(13); clf; 
loglog(freq,smooth(f1,sf),'DisplayName','Left');hold on;
f2=(abs(fft(d1.dd.data{2})))./sqrt(size(d1.dd.data{2},1)*nyquist*2); f2=f2(1:end/2); f2(1)=0;
loglog(freq,smooth(f2,sf),'r','DisplayName','Right');
legend show
title('Noise on 2013-05-18'); xlabel('Frequency (Hz)'); ylabel('Noise $$(V/\sqrt{Hz})$$','Interpreter','latex','FontSize',13);




