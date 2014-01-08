%% 
files = get_files('sm*.mat');
%% this cell and next for gain sweep with true rotating frame ramsey
data = []; T2= []; yv = [];
for j = 1:length(files)
    d = ana_avg(files{j});
    yv(j) = d.scan.data.CSGain;
    data(j,:) = squeeze(mean(d.data{1}(:,:,1001:end)));
    pars = fitoscillations(d.xv{1}(1001:end),data(j,:),'gauss','plinit plfit');
    T2(j) = pars(5);
end
xv = d.xv{1}(1001:end);

%%
close all
figure(1); clf;
imagesc(xv,yv,data);
xlabel('evolution time');
ylabel('CDS gain');
title('rotating ramsey with CDS');
figure(2); clf; hold on;
c = 'rgbcmyk'; c = [c c c c c c c c c c c];
legstr = {};
for j = 1:2:length(T2)
    plot(xv,data(j,:),c(j),'DisplayName',sprintf('DAC gain = %d',yv(j)));
end
xlabel('evolution time (ns)');
ylabel('Triplet Probability');

legend show
figure(3); clf;
plot(yv,T2,'x','MarkerSize',10,'LineWidth',4);
xlabel('CDS gain');
ylabel('T2* (ns)');


%% this cell and next for offset sweep with true rotating frame ramsey
data = []; T2= []; yv = [];
for j = 1:length(files)
    d = ana_avg(files{j});
    blah=sscanf(d.scan.data.prettyname,'dbz_rotramsey_VCO_%f_R_%i'); yv(j) = blah(1);
    %yv(j) = d.scan.data.CSOffset;
    data(j,:) = squeeze(mean(d.data{1}));
    %pars = fitoscillations(d.xv{1}(1001:end),data(j,:),'gauss','plinit plfit');
    %T2(j) = pars(5);
end
xv = d.xv{1};
%%
figure(1); clf; imagesc(xv,yv,data);
xlabel('evolution time (ns)');
ylabel('VCO control voltage (V)');
title('offfset scan');
figure(2); clf; hold on
c = 'rgbcmyk'; c = [c c c c];
for j = 1:length(files)
    plot(xv,data(j,:),c(j),'DisplayName',sprintf('Offset = %.2fV',yv(j)));
end
xlabel('evolution time');
ylabel('Triplet Probability');
legend show