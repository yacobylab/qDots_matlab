function [mvis mtime vthresh t1 Tfit Shist Thist Vt] = anaRead2(scan, data)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


pulselength = plsinfo('zl', scan.data.pulsegroups(1).name);
pulselength = abs(pulselength(1));

dt = 1/scan.configfn.args{3}(2);

data1 = data{1};
data2 = reshape(data1',1e-9*2*pulselength/dt,size(data1,1)*size(data1,2)/(1e-9*2*pulselength/dt));

Smask = 1:120;
Tmask = 1:120;

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
params = fitwrap('',dt*(1:length(avgT)), avgT', beta0, fitfcn,[ 0 1 1]);
t1=params(3);
figure(80); clf; hold on;
subplot(4,2,7); hold on;
plot(dt*(1:length(avgT)), avgT',dt*(1:length(avgT)), fitfcn(params, dt*(1:length(avgT))));
title(sprintf('T_{1} = %.2f \\mus', 1e6*t1));


singHist = cumsum(singlet,1)./repmat((1:size(singlet,1))',1,size(singlet,2)); %poor man's average
tripHist = cumsum(triplet,1)./repmat((1:size(triplet,1))',1,size(triplet,2)); 


allData = [singHist tripHist];
cen = mean(allData(:));
rng = 3*std(allData(:));
Vt=linspace(cen-rng,cen+rng,512);

Shist = histc(singHist', Vt);
subplot(2,2,1); hold on;
imagesc(1e6*T, Vt,Shist); title('Singlet Histogram');
xlabel('T_{meas}(\mus)');

Thist = histc(tripHist', Vt);
subplot(2,2,2); hold on;
imagesc(1e6*T, Vt,Thist); title('Triplet Histogram');
xlabel('T_{meas}(\mus)');

% Try to fit some histograms
for i=1:size(Thist,2)    
  [fitfn, initfn] = getfn(T(i)/t1);
  Tfit(i,:)=fitwrap('ause linit plfit', [Vt  Vt], [Shist(:,i) ; Thist(:,i)]', initfn, fitfn, [1 1 1 1 1 1 0 0 1]);
end
%
sfid = 1-cumsum(Shist,1)./repmat(sum(Shist,1),length(Vt),1);
tfid = 1-cumsum(Thist(end:-1:1,:),1)./repmat(sum(Thist,1),length(Vt),1);
tfid = tfid(end:-1:1,:);
[maxsfid maxsfid_i] = max(sfid);
[maxtfid maxtfid_i] = max(tfid);
vis = 1-tfid-sfid;
[maxvis maxvis_i] = max(vis);

[visibility Tmeas] = max(maxvis);

subplot(2,2,4); hold on;
plot(vis);
title(sprintf('Visibility = %.2f  T_{meas} = %.2f \\musec\n V_{threshold} = %.1f mV ', visibility, 1e6*dt*Tmeas,1e3*Vt(maxvis_i(Tmeas))));

subplot(4,2,5);
plot(1e6*T,maxvis);
xlabel('T (\mus)');
title('Best Visibility');

fprintf('Visibility = %.2f (1-%g)\n',visibility,1-visibility)
fprintf('Tmeas = %.2d usec\n Vthrehold = %3d \n',dt*Tmeas*1e6,Vt(maxvis_i(Tmeas)));
fprintf('T1 = %g us\n',1e6*params(3));
mvis=visibility;
mtime = dt*Tmeas;
vthresh=Vt(maxvis_i(Tmeas));
t1=1e6*params(3);
return;

function [cfitfn, cinitfn] = getfn(t1)

distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
     a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
     (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [t_meas/T1, rms amp noise/peak spacing]
 
fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [ center between peaks, 1/spacing, coeff left peak, coeff right peak, t_m/T1,left, t_m/T1,right, rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

% parameters [centering, 1/spacing 1s 1t 2s 2t tm/t1_l tm/t1_r sig/spacing]
cfitfn = @(a, x) [fitfn(a([1 2 3 4 7 8 9]),x(1:end/2)) fitfn(a([1 2 5 6 7 8 9]),x((end/2+1):end))];
cinitfn.fn = @(x, y)[sum(x.*y)/sum(y), 125, max(y), max(y), max(y), max(y), 1e-10, t1, .2];
cinitfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];

return;
