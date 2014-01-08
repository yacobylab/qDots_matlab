function [mvis mtime vthresh t1 Shist Thist T Vt] = anaRead1(scan, data,fignum)
%function [mvis mtime vthresh t1] = anaRead1(scan, data,fignum)
% returs maximum visibility, optimal measurement time, optimal Vthreshold
% and computes T1 a la C. Barthel paper on singleshot readout: Vis =
% Sfid+Tfid-1.
% Takes scan and data (scan loads singlet and measures (oversampled) and
% then triplet and measures (oversampled).  
% format of data is [S T S T....] (100 of each per row). for now it is 50
% of such rows, but code should handle different sizes. 
% Shist and Thist are all the singlet and triplet histograms
% T is the vector of times for Thist and Shist
% Vt is the vector of RF voltages.
% fignum is the figure number
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if ~exist('fignum','var')
    fignum=79;
end
pulselength = plsinfo('zl', scan.data.pulsegroups(1).name);
pulselength = abs(pulselength(1));

dt = 1/scan.configfn.args{3}(2); %integration time bin

data1 = data{1};
data2 = reshape(data1',1e-9*2*pulselength/dt,size(data1,1)*size(data1,2)/(1e-9*2*pulselength/dt));

% New; make mask robust against changes in sampling rate.
%mask for getting rid of manipulation time. hardcoded for now
drop_time=120*20e-9;
Smask = 1:ceil(drop_time/dt); 
Tmask = 1:ceil(drop_time/dt);
max(Smask);

singlet = data2(1:size(data2,1)/2,:);
triplet = data2(1+size(data2,1)/2:size(data2,1),:);
singlet(Smask,:)= [];
triplet(Tmask,:)= [];
T = dt*(1:size(singlet,1));

%calc T1
%avgS = mean(singlet,2);
avgT = mean(triplet,2)- mean(singlet,2);
beta0 = [0, range(avgT), 1/abs(range(avgT))*dt];
fitfcn = @(p,x)p(1)+p(2).*exp(-1*x./p(3));
params = fitwrap('',dt*(1:length(avgT)), avgT', beta0, fitfcn, [ 0 1 1]);
figure(fignum); clf; hold on;
subplot(4,2,7); hold on;
plot(dt*(1:length(avgT)), avgT',dt*(1:length(avgT)), fitfcn(params, dt*(1:length(avgT))));
title(sprintf('T_{1} = %.2f \\mus', 1e6*params(3)));


singHist = cumsum(singlet,1)./repmat((1:size(singlet,1))',1,size(singlet,2)); %poor man's average
tripHist = cumsum(triplet,1)./repmat((1:size(triplet,1))',1,size(triplet,2)); 
allData = [singHist tripHist];
cen = mean(allData(:));
rng = 3*std(allData(:));
Vt=linspace(cen-rng,cen+rng,512); %512 voltage bins, hardcoded for now.


Shist = histc(singHist', Vt);
subplot(2,2,1); hold on;
imagesc(1e6*T, Vt,log(Shist+1)); title('Singlet Histogram');
xlabel('T_{meas}(\mus)');

Thist = histc(tripHist', Vt);
subplot(2,2,2); hold on;
imagesc(1e6*T, Vt,log(Thist+1)); title('Triplet Histogram');
xlabel('T_{meas}(\mus)');

%fid = 1-(#misidentified/#tot); cumsum is good way of getting this.
sfid = 1-cumsum(Shist,1)./repmat(sum(Shist,1),length(Vt),1);
tfid = 1-cumsum(Thist(end:-1:1,:),1)./repmat(sum(Thist,1),length(Vt),1);
tfid = tfid(end:-1:1,:);
vis = 1-tfid-sfid;
[maxvis maxvis_i] = max(vis); %maxvis_i is optimal Vthreshold

[visibility Tmeas] = max(maxvis)

subplot(2,2,4); hold on;
plot(vis);
title(sprintf('Visibility = %.2f  T_{meas} = %.2f \\musec\n V_{threshold} = %.1f mV ', visibility, 1e6*dt*Tmeas,1e3*Vt(maxvis_i(Tmeas))));

subplot(4,2,5);
plot(1e6*T,maxvis);
xlabel('T (\mus)');
title('Best Visibility');

fprintf('Visibility = %.2f \n Tmeas = %.2d usec\n Vthreshold = %3d \n', visibility, dt*Tmeas*1e6,Vt(maxvis_i(Tmeas)));
fprintf('V_s = %.2f, V_t = %.2f\n',1-sfid(maxvis_i(Tmeas),Tmeas),1-tfid(maxvis_i(Tmeas),Tmeas));
fprintf('T1 = %g us\n',1e6*params(3));
mvis=visibility;
mtime = dt*Tmeas;
vthresh=Vt(maxvis_i(Tmeas));
t1=1e6*params(3);
return;
