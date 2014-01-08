file = uigetfile('sm*.mat','MultiSelect','on');
file = sort(file);
%     'sm_constrast_L_2087.mat'
%     'sm_constrast_L_2088.mat'
%     'sm_constrast_L_2089.mat'
%     'sm_constrast_L_2090.mat'
%     'sm_constrast_L_2091.mat'
%%
s =ana_avg(file,struct('opts','hd'));
figure(1); 
legstr = {'dBz', 'echo (dbz)', 'echo (adprep)', 'ramsey (dbz)', 'ramsey (adprep)'};
legend(legstr);
close all
%% plot dbzs
c = 'rgbcmyk';
c = [c c c c c];
legstr2 = {'No HW, 1ns SW', 'No HW No SW', '300MHz HW, no SW',...
    '300MHz HW, 1ns SW','300MHz HW, 3ns SW', 'No HW no SW','Compensated'};
for j = 1:length(file)
   legstr3{j}= [legstr2{j}, ' ' file{j}]; 
end
ind = 1;
figure(2); clf; hold on;
for j = 1:length(s)
   plot(squeeze(mean(s(j).data{1}(:,ind,:))),c(j));    
end

legend(legstr3);
title('dBz');

%% RamseyE (dbz)

ind = 2;
figure(3); clf; hold on;
for j = 1:length(s)
   plot(squeeze(mean(s(j).data{1}(:,ind,:))),c(j));    
end

legend(legstr3);
title('echo (dbz)');

%% RamseyE (adprep)

ind = 3;
figure(4); clf; hold on;
for j = 1:length(s)
   plot(squeeze(mean(s(j).data{1}(:,ind,:))),c(j));    
end

legend(legstr3);
title('echo (adprep)');