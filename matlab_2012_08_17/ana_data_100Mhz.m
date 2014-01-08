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