%% 
clear dd;
smset('PulseLine',awgseqind('chrg_1_LR'));
dummy = sm_getDAQBuffers([1 2],100e6,2e6); %20ms of data
dd.data = dummy;
dd.time =now;
dd.note = 'grounded';
fname = smnext('sm_100MHz_LR');
save(fname,'dd');

%%
n=50
%for i=1:50 
 
%% 2013_05_18
clear dd;
smset('PulseLine',awgseqind('chrg_1_LR'));
dummy = sm_getDAQBuffers([1 2],100e6,2e6); %20ms of data
dd.data = dummy;
dd.time =now;
dd.note = 'sensitive';
fname = smnext('sm_100MHz_LR');
save(fname,'dd');
