%%
load('results_2013_07_16_SBAD_no_ana.mat');
%%
ftime = 10:10:200;
mdbz = .089;
results = results2;
for j = 1:length(results2)
    dbz = get_dbz(results2(j).fourier_data,ftime, mdbz);
    results(j).dbzs=dbz;
    ti = linspace(.95*results(j).t(1),results(j).t(end)*1.05,length(results(j).t));
    for k = 1:size(results(j).evo_data,1)
        ts = results(j).t*results(j).dbzs(k)/mean(results(j).dbzs);
       datatmp(k,:) = interp1(ts,results(j).evo_data(k,:),ti);
    end
    results(j).newdata = nanmean(datatmp);
    results(j).newtime = ti;
end

%% plot the data

figure(3); clf; hold on;
for j = 1:length(results)
    mask2 = 1:length(results(1).newtime);
    plot(results(j).newtime(mask2),results(j).newdata(mask2));    
end
%%
%%
%% for sydney team
clear results dataI
alphabest = 0;%d_good.alphabest;
four_inds = 1:200;
evo_inds = (1+four_inds(end)):250;
for j = 1:length(files)
    %d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt flip','mdbz',.089,'shot_include',1:75));
    d2 = ana_rescale_dbz_v3(files{j},struct('opts','quick noppt flip','mdbz',.089));
    results(j).dbzs = d2.pmax;
    d=ana_avg(files{j},'noppt'); close(2);
    results(j).file = d.filename;
    results(j).fourier_data = squeeze(d.data{1}(:,four_inds));
    results(j).evo_data = squeeze(d.data{1}(:,evo_inds));
    data = squeeze(d.data{1});
    t = d.xv{1}(evo_inds);
    results(j).t = t;
    ti = linspace(.95*t(1),t(end)*1.05,length(t));
    %ph_good = linspace(min(t(1)*results(j).pmax)-pi,pi+max(t(end)*results(j).pmax),length(t));
    %factor = 5;
    for k = 1:size(data,1)        
        %tnew = linspace(ts(1),ts(end),length(ts)*factor);
        %ph = results(j).scale(k)*t;
        %dataI(k,:)=interp1(ph,data(k,evo_inds),ph_good); 
        %ts = results(j).tfunc(t,0,results(j).scale(k));
        ts = t*results(j).dbzs(k)/mean(results(j).dbzs);
        if all(ts==0)
           ts = ti; 
        end
        dataI(k,:)=interp1(ts,data(k,evo_inds),ti); 
    end
    mmask = 1:size(dataI,2);
    results(j).newdata = nanmean(dataI(:,mmask));
    results(j).newtime = ti(mmask);
    %keyboard
    fprintf('done with %i of %i \n',j,length(files));
end
close all;