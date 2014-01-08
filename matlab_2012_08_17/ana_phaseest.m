s=ana_avg;
%%
data = squeeze(s.data{1});
svals = mean(data(:,1:end/2),2);
meanS = mean(svals);
svals = svals-meanS;
dbzval = data(:,1+end/2:end);
meandbz = mean(dbzval);
%%
figure(1); clf; 
plot(meandbz)
figure(2); clf; 
plot(svals,'.');
figure(3); clf;
imagesc(dbzval)
%% fit dbz avg
tv = (1:length(meandbz))*1e-9;
fitfn = @(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2);
beta0 = [mean(meandbz) range(meandbz), .001,  2*pi*75e6,  1e-9, 100e5];
pars = fitwrap('plinit plfit',tv, meandbz, beta0, fitfn);

%%
mdbz = pars(4)/(2*pi);
tevo = 17e-9;
ddbzs = -1*meanS./svals;
frqs = 2*pi*mdbz + svals.*sqrt((1-meanS^2))/(tevo);
newtv = repmat(tv,size(data,1),1);
for j = 1:length(newtv)
    newtv(j,:) = newtv(j,:)*2*pi*(mdbz/frqs(j));
end

figure(4); clf; hold on;
for j = 1:size(dbzval,1)
   plot(newtv(j,:),dbzval(j,:),'.'); 
end


%%
s=ana_avg; close all;
data=squeeze(s.data{1});
d1= mean(data(:,1:100),2);
d2= mean(data(:,101:200),2);
close all
plot(d1,d2,'.')