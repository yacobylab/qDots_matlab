%% conditional evolution with control at different epsilons

files=(uigetfile('sm*.mat','MultiSelect','on')); 
if ~iscell(files)
   files = {files}; 
end
% sm_condEvow_LR_1682.mat 
% sm_condEvow_LR_1681.mat 
% sm_condEvow_LR_1680.mat 
% sm_condEvow_LR_1679.mat 
% sm_condEvow_LR_1678.mat 
% sm_condEvow_LR_1677.mat 
% sm_condEvow_LR_1676.mat 

%%
tv = (1:65);
c ='rgbcmyk'; c = [c c c c c c];
results = [];
for j = 1:length(files)
    s= load(files{j});
    scantime =getscantime(s.scan,s.data);
    plsgrps = plsinfo('gd',s.scan.data.pulsegroups(2).name);
    params=plsinfo('params', plsgrps.pulses.groups{2},[],scantime);
    
    dataL = squeeze(nanmean(s.data{1}));
    dataR = squeeze(nanmean(s.data{2}));
    dbz = dataL(1,:);
    dataL = dataL(2:end,:); % dont want the dbz;
    dataR = dataR(2:end,:); % dont want the dbz;
    
    offset = 0;
    bestJ = [];
    for k = 2:size(dataR,1)
       bestJ(k) = 1e3*guessfreq(tv,dataR(k,:))/(2*pi);
    end
    
    figure(2); clf; hold on;
    for k = 1:size(dataL,1)
        plot(tv,dataL(k,:)+offset, c(k));
        offset = offset+1e-3;
    end
    
    xlabel('\taux_J');
    ylabel('\tau_{dBz}');
    title('Conditional Evolution');
    
    tauJ = ginput(1);
    results(end+1).tauEnt = round(tauJ(1));
    results(end).eps = params(2);
    results(end).bestJ = bestJ;
    
end
figure(3); clf; 
plot([results.eps],[results.tauEnt],'x','MarkerSize',8);
xlabel('\epsilon_{control} (mV)');
ylabel('first beat time (ns)');
%% same as above but for a left control. 

figure(4); clf;
for j =1:length(results)
 plot(results(j).bestJ(2:end),'.')
 pause
end


tv = (1:65);
c ='rgbcmyk'; c = [c c c c c c];
results = [];
for j = 1:length(files)
    s= load(files{j});
    scantime =getscantime(s.scan,s.data);
    plsgrps = plsinfo('gd',s.scan.data.pulsegroups(2).name);
    params=plsinfo('params', plsgrps.pulses.groups{2},[],scantime);
    
    dataL = squeeze(nanmean(s.data{1}));
    dataR = squeeze(nanmean(s.data{2}));
    dbz = dataL(1,:);
    dataL = dataL(2:end,:); % dont want the dbz;
    dataR = dataR(2:end,:); % dont want the dbz;
    
    offset = 0;
    bestJ = [];
    for k = 2:size(dataL,1)
       bestJ(k) = 1e3*guessfreq(tv,dataL(k,:))/(2*pi);
    end
    
    figure(2); clf; hold on;
    for k = 1:size(dataL,1)
        plot(tv,dataR(k,:)+offset, c(k));
        offset = offset+1e-3;
    end
    
    xlabel('\taux_J');
    ylabel('\tau_{dBz}');
    title('Conditional Evolution');
    
    tauJ = ginput(1);
    results(end+1).tauEnt = round(tauJ(1));
    results(end).eps = params(2);
    results(end).bestJ = bestJ;
    
end
figure(3); clf; 
plot([results.eps],[results.tauEnt],'x','MarkerSize',8);
xlabel('\epsilon_{control} (mV)');
ylabel('first beat time (ns)');