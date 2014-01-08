%% for dbz ramsey correction 
%
d = ana_avg;
gain=[];
close all;
figure(666); clf; hold on; c = 'krgbcmy';
results = struct();
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^2);
for j =1:length(d)
    cnum=mod(j,7)+1;
    results(j).data = squeeze(nanmean(d(j).data{1}(:,:,201:end)))';
    results(j).xv = d(j).xv{1}(201:end);
    figure(666);
    plot(results(j).xv,results(j).data,c(cnum),'DisplayName',d(j).filename(4:end-4));    
    beta0 = [0.4,0.3,200,.1,700];
    gain(j)=d(j).scan.data.FPGA.Gain; 
    %gain(j)=d(j).scan.data.FPGA.Offset
    try
       pars = fitwrap('plinit plfit',results(j).xv,results(j).data,beta0,fitfn);
       results(j).T2 =abs(pars(5));
       results(j).fitpars = pars;       
    catch
        results(j).T2 = NaN;
        results(j).fitpars = [];
    end
end
figure(666);
hleg=legend('toggle');
set(hleg,'Interpreter','none');


[~,ind]=max([results.T2]);
%figure(667); clf; plot([results.T2],'x','LineWidth',5,'MarkerSize',5)
figure(667); clf; plot(gain,[results.T2],'x','LineWidth',5,'MarkerSize',5)
xlabel('FPGA Gain'); 
%xlabel('FPGA Offset')
ylabel('T2*') 
%%
figure(2); clf; hold on
plot(results(ind).xv,results(ind).data,'b.');
plot(results(ind).xv,fitfn(results(ind).fitpars,results(ind).xv),'r')
xlabel('evolution time (ns)');
ylabel('Ramsey amplitude');
title(sprintf('T_2{}^* = %.1fns',results(ind).T2))

%%
d= ana_avg; close all
results2 = struct;
results2.data = [];
results2.xv = [];
for j = 1:length(d)
    results2.xv = [results2.xv, d(j).xv{1}(201:end)];
    results2.data = [results2.data, squeeze(nanmean(d(j).data{1}(:,:,201:end)))'];
end
figure(661); clf;
plot(results2.xv,results2.data,'k.')
%% For ramsey ramsey correction. 
d = ana_avg;
gain=[];
close all;
figure(666); clf; hold on; c = 'krgbcmy';
results = struct();
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^p(6));
clear T2err; clear T2rng; clear results;  
for j =1:length(d)
   gd=plsinfo('gd',d(j).scan.data.pulsegroups(1).name,[],d(j).scantime);
    fourier_size=size(gd.varpar{1},1);
    evo_size=size(gd.varpar{2},1);
    evoinds=fourier_size+[1:evo_size]; 
    cnum=mod(j,7)+1;
    results(j).data = squeeze(nanmean(d(j).data{1}(:,:,evoinds)))';
    results(j).xv = d(j).xv{1}(evoinds);
    plot(results(j).xv,results(j).data,c(cnum),'DisplayName',d(j).filename(4:end-4));    
    pars= [0.4,0.2,110,.1,700,1.5];
    if isfield(d(j).scan.data.FPGA,'Gain') && ~isempty(d(j).scan.data.FPGA.Gain)
        gain(j)=d(j).scan.data.FPGA.Gain;
    else
        gain(j)=nan;
    end
    %gain(j)=d(j).scan.data.FPGA.Offset
    try
      pars = fitwrap('plinit plfit',results(j).xv,results(j).data,pars,fitfn,[1 1 1 1 1 0]);
       results(j).pars=pars; 
      [pars,res,jac,sig] = fitwrap('plinit plfit',results(j).xv,results(j).data,pars,fitfn,[1 1 1 1 1 1]);
      jac=jac';
      jac=reshape(jac,length(jac)/6,6); 
        
      ci=nlparci(pars,res,'jacobian',jac);
       % ci=nlparci(pars,res,'covar',sig)
       results(j).T2 =abs(pars(5));
       results(j).T2rng=ci(5,:)'; 
       results(j).fitpars = pars;       
       results(j).exp=pars(6); 
       results(j).exprng=ci(6,:); 
    catch
        results(j).T2 = NaN;
        results(j).fitpars = [];
    end
end

figure(666);
hleg=legend('toggle');
set(hleg,'Interpreter','none');
[~,ind]=max([results.T2]);
% figure(667); clf; plot(gain,[results.T2],'x','LineWidth',3,'MarkerSize',5); hold on; 
T2rng=[results.T2rng];
% plot(gain,T2rng(1,:),gain,T2rng(2,:),'x-','MarkerSize',5)
figure(667); clf;
T2err(1,:)=[results.T2]-T2rng(1,:); 
T2err(2,:)=-([results.T2]-T2rng(2,:)); 
errorbar(gain,[results.T2],T2err(1,:),T2err(2,:),'x-','LineWidth',1,'MarkerSize',10); 
xlabel('FPGA Gain'); 
%xlabel('FPGA Offset')
ylabel('T2*') 


figure(2); clf; hold on
plot(results(ind).xv,results(ind).data,'b.-');
plot(results(ind).xv,fitfn(results(ind).fitpars,results(ind).xv),'r')
xlabel('evolution time (ns)');
ylabel('Ramsey amplitude');
title(sprintf('T_2{}^* = %1.1fns, exp=%1.1f (%1.1f - %1.1f)',results(ind).T2,results(ind).exp,results(ind).exprng(1),results(ind).exprng(2)))


fpga_freqs=d(ind).scan.data.FPGA.freqs; 
fpga_data=d(ind).scan.data.FPGA.data; 
figure(11); 
[ct,bins]=hist(fpga_freqs,50); 
plot(bins,ct);
xlabel('Frequency Index') 
ylabel('Count')

%% hack something together for multiple pulse lengths. 


d = ana_avg;
gain=[];
%close all;
figure(666); clf; hold on; c = 'krgbcmy';
results = struct();
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^p(6));
clear T2err; clear T2rng; clear results;  
for j =1:length(d)
    gd=plsinfo('gd',d(j).scan.data.pulsegroups(1).name,[],d(j).scantime);
    evo_size=0; 
    evo_st=find(gd.pulses==72);
    evo_st=evo_st(1); 
    varlengths=cellfun(@length,gd.varpar); 
    fourier_size=sum(varlengths(1:evo_st-1)); 
    evo_size=sum(varlengths(evo_st:end)); 
    evoinds=fourier_size+[6:evo_size]; 
    cnum=mod(j,7)+1;
    results(j).data = squeeze(nanmean(d(j).data{1}(:,:,evoinds)))';
    results(j).xv = d(j).xv{1}(evoinds);
    plot(results(j).xv,results(j).data,c(cnum),'DisplayName',d(j).filename(4:end-4));    
    pars= [0.4,0.2,110,.1,700,1.5];
    gain(j)=d(j).scan.data.FPGA.Gain; 
    %gain(j)=d(j).scan.data.FPGA.Offset
    try
      pars = fitwrap('plinit plfit',results(j).xv,results(j).data,pars,fitfn,[1 1 1 1 1 0]);
      [pars,res,jac,sig] = fitwrap('plinit plfit',results(j).xv,results(j).data,pars,fitfn,[1 1 1 1 1 1]);
      jac=jac';
      jac=reshape(jac,length(jac)/6,6); 
        
      ci=nlparci(pars,res,'jacobian',jac);
       % ci=nlparci(pars,res,'covar',sig)
       results(j).T2 =abs(pars(5));
       results(j).T2rng=ci(5,:)'; 
       results(j).fitpars = pars;       
       results(j).exp=pars(6); 
       results(j).exprng=ci(6,:); 
    catch
        results(j).T2 = NaN;
        results(j).fitpars = [];
    end
end


%% plot evos split into multiple data sets. 
d = ana_avg;
d=flipud(d)
gain=[];
%close all;
figure(668); clf; hold on; 
figure(666); clf; hold on; c = 'krgbcmy';
results = struct();
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^p(6));
clear T2err; clear T2rng; clear results;  
c='bbbbbbb';
ydata=[]; xdata=[];
for j =1:length(d)
    gd=plsinfo('gd',d(j).scan.data.pulsegroups(1).name,[],d(j).scantime);
    evo_size=0; 
    evo_st=find(gd.pulses==72); 
    varlengths=cellfun(@length,gd.varpar); 
    fourier_size=sum(varlengths(1:evo_st-1)); 
    evo_size=sum(varlengths(evo_st:end)); 
    evoinds=fourier_size+[6:evo_size]; 
    cnum=mod(j,7)+1;
    results(j).data = squeeze(nanmean(d(j).data{1}(:,:,evoinds)))';
    results(j).xv = d(j).xv{1}(evoinds);
    figure(666);
    h=plot(results(j).xv,results(j).data,'.-','Color',c(cnum),'DisplayName',d(j).filename(4:end-4));        
%     fpga_freqs=d(j).scan.data.FPGA.freqs;
%     figure(668);
%     inds=find(fpga_freqs>128 | fpga_freqs<=3);
%     datafilt=d(j).data{1};
%     datafilt(inds,:,:)=[];
%     datafiltm=squeeze(nanmean(datafilt(:,:,evoinds)));
%     plot(results(j).xv,datafiltm,'.-')
    beta0= [0.4,0.2,110,-.1,1200,1.5];
    %gain(j)=d(j).scan.data.FPGA.Gain; 
    %assumes you go in revers order. 
    ydata=[ydata, fliplr(results(j).data)]; 
    xdata=[xdata, fliplr(results(j).xv)]; 
end

ydata=fliplr(ydata); 
xdata=fliplr(xdata); 
%%
fpga_freqs=d(1).scan.data.FPGA.freqs; 
figure(11); clf;
[ct,bins]=hist(fpga_freqs,50); 
plot(bins,ct);
xlabel('Frequency Index') 
ylabel('Count')

inds=find(fpga_freqs>128)
results(j).data(inds,:)=[]
figure 668; clf; 

%% combine all data together and perform one fit. 
fitfn = @(p,x) p(1)+p(2)*cos(2*pi*x/p(3)+p(4)).*exp(-(x/p(5)).^p(6));
     try
      pars = fitwrap('plinit plfit',xdata,ydata,beta0,fitfn,[1 1 1 1 1 1]);
      [pars,res,jac,sig] = fitwrap('plfit',xdata,ydata,pars,fitfn,[1 1 1 1 1 1]);
      jac=jac';
      jac=reshape(jac,length(jac)/6,6);         
      ci=nlparci(pars,res,'jacobian',jac);
      T2 =abs(pars(5));
      T2rng=ci(5,:)'; 
      fitpars = pars;       
      expnt=pars(6); 
      exprng=ci(6,:); 
    catch
        results(j).T2 = NaN;
        results(j).fitpars = [];
    end

title(sprintf('T2* is %.0f ns and decay exp is %0.1f (%0.1f to %0.1f)',T2,expnt,exprng(1),exprng(2))); 

%% perform an mfit. 
data = [];
model = [];
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(7))).*exp(-((x+p(5))/p(4)).^p(6));
mask = [1 1 1 1 0 0];
beta0 = [.36, .36, 150, 500, 0, 1.5];
boundary = 0;
for j = 1:length(results)
    if all((~isnan(results(j).data)))
        maskfit=6:length(results(j).xv);
        %maskfit=1:length(results(j).xv); 
        data(end+1).y = results(j).data(:,maskfit); 
        data(end).x = results(j).xv(:,maskfit);
        model(end+1).fn = fitfn;
        model(end).pt = @(p) [p(1:6), p(7+boundary)];
        boundary = boundary+1;
        mask = [mask, 1];
        beta0 = [beta0, .01];
    end
end

pars = mfitwrap(data,model,beta0,'plinit plfit lm',mask);
mask(6) = 1;
pars = mfitwrap(data,model,pars,'plinit plfit lm',mask);
figure(666)
    for j=1:length(results); 
        %plot(data(j).x,out.results(j).data(ind,mask2),'x-'); hold on; 
        plot(data(j).x,fitfn(model(j).pt(pars),data(j).x),'r');
    end
title(sprintf('T2* is %1.1f with exponent %1.1f',pars(4),pars(6))); 




