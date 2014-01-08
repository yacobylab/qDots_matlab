sQ = [];
for j = 1:64
PauliPlot(dataexp(j,:));
pause

sQ(j) = sum(dataexp(j,1:6).^2);

end

%%
figure(2);
clf;
plot(sQ)

%%
blah =ana_TomoUnpack;
sQ2= []; tQ2 = [];
c = 'rgbcmyk';
figure(3); clf; hold on;

for j = 1:size(blah,2)
   PauliPlot(blah(:,j));
   sQ2(j) = sum(blah(1:6,j).^2);
   tQ2(j) = sum(blah(7:end,j).^2);
end
dt = .01;
tvals = dt:dt:64*dt;
mask = 3:length(tvals);
figure(3); %clf;
plot(tvals(mask),sQ2(mask),'-.'); hold on;
plot(tvals(mask),tQ2(mask),'rx-');
legend('Single Qubit Weight','Two qubit weight');
%%
figure(989); clf; hold on;
bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};
colors=[1 0 0; 0 0 1 ; .7 0 .7];
cind=[1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 ];
l=15;
for i=1:l
    %area(fitdata2.mfdata(1).x, fitdata2.dataexp(i,:)+2*(l-i), 2*(l-i),'FaceColor',colors(cind(i),:));
    area(tvals, blah(i,:)+2*(l-i), 2*(l-i),'FaceColor',colors(cind(i),:));
end
set(gca,'YTick',2*(0:1:(l-1)));
set(gca,'YTickLabel',bs(end:-1:1));
set(gca,'YLim',[-1 l*2-1]);
xlabel('T_{entangle} (\mus)');
%%
figure(3); clf;
plot(1e3*tvals,blah([2 9],:)); 

%%
ss= load(uigetfile('sm*.mat'));
ddL = squeeze(nanmean(s.data{1}));
ddR = squeeze(nanmean(s.data{2}));
figure(9); clf;
imagesc(ddL);
figure(10); clf; 
imagesc(ddR)