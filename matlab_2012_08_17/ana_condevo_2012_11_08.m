%%
s=load(uigetfile('sm*.mat')); % sm_cond_evo_LR_1453.mat
dataL = squeeze(nanmean(s.data{1}));
dataR = squeeze(nanmean(s.data{2}));
dbz = dataL(1,:);
dataL = dataL(2:end,:); % dont want the dbz;
dataR = dataR(2:end,:); % dont want the dbz;


%%
tv = (1:65);
offset = 0;
c ='rgbcmyk'; c = [c c c c c c];
figure(2); clf; hold on;
for j = 1:size(dataL,1)
    plot(tv,dataR(j,:)+offset, c(j));
    offset = offset+1e-3;
end

xlabel('\tau_J');
ylabel('\tau_{dBz}');
title('Conditional Evolution');

%% anaylze raw mode data

ss= load(uigetfile('sm*.mat')); %sm_cond_evo_LR_1454.mat
rdataL = ss.data{4};
rdataR = ss.data{5};
% unpacked as 50 x 34 x 3250 elements
% = nrep, plsgrp, nloop*(number of pulses in group)

%% This cell makes you decide which pulse group (dbz_evo time) and which J_evo time to use

ind = 7; % which pulsegroup. remember, there is a dbz reference > ind = 1
tauJ = 15; % the J-evo time to look at

lcut = squeeze(rdataL(:,ind,:));
lcut = reshape(lcut,ss.scan.data.conf.nrep,ss.scan.data.conf.nloop,size(ss.data{1},3));

hhL = (lcut(:,:,tauJ));

rcut = squeeze(rdataR(:,ind,:));
rcut = reshape(rcut,ss.scan.data.conf.nrep,ss.scan.data.conf.nloop,size(ss.data{1},3));

hhR = rcut(:,:,tauJ);

figure(2); clf; hold on;
plot(hhL,hhR,'.');
xlabel('Left Sensor (target qubit)');
ylabel('Right Sensor (control qubit)');

%% sweep through all of the J_evo times (you pick dbz time

ind = 9; % which pulsegroup. remember, there is a dbz reference > ind = 1

lcut = squeeze(rdataL(:,ind,:));
lcut = reshape(lcut,ss.scan.data.conf.nrep,ss.scan.data.conf.nloop,size(ss.data{1},3));

rcut = squeeze(rdataR(:,ind,:));
rcut = reshape(rcut,ss.scan.data.conf.nrep,ss.scan.data.conf.nloop,size(ss.data{1},3));

figure(3); clf;

for j = 1:size(lcut,3)
    plot(lcut(:,:,j),rcut(:,:,j),'.')
    xlabel('Left Sensor (target qubit)');
    ylabel('Right Sensor (control qubit)');
    pause
end


%% conditional evolution with control at different epsilons

s=load(uigetfile('sm*.mat')); % sm_cond_evo_LR_1453.mat
dataL = squeeze(nanmean(s.data{1}));
dataR = squeeze(nanmean(s.data{2}));
dbz = dataL(1,:);
dataL = dataL(2:end,:); % dont want the dbz;
dataR = dataR(2:end,:); % dont want the dbz;



tv = (1:65);
offset = 0;
c ='rgbcmyk'; c = [c c c c c c];
figure(2); clf; hold on;
for j = 1:size(dataL,1)
    plot(tv,dataR(j,:)+offset, c(j));
    offset = offset+1e-3;
end

xlabel('\tau_J');
ylabel('\tau_{dBz}');
title('Conditional Evolution');

%% 
eps_cntl = [0 .5 1 1.4 1.8];
tau_beat =[16 14 14 15 14];
