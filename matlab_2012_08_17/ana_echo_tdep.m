%% get the files
%files={'sm_RamseyE_5k_R_1999.mat','sm_RamseyE_R_5k_2003.mat','sm_RamseyE_R_5k_2004.mat','sm_RamseyE_R_5k_2005.mat','sm_RmseyE_R_5k_2006.mat','sm_RamseyE_R_5k_2007.mat','sm_RamseyE_R_5k_2008.mat','sm_RamseyE_R_5k_2009.mat','sm_RamseyE_R_5k_2010.mat','sm_RamseyE_R_5k_2011.mat','sm_RamseyE_R_5k_2012.mat','sm_RamseyE_R_5k_2013.mat','sm_RamseyE_R_5k_2014.mat','sm_RamseyE_R_5k_2015.mat','sm_RamseyE_R_5k_2016.mat','sm_RamseyE_R_5k_2017.mat','sm_RamseyE_R_5k_2018.mat','sm_RamseyE_R_5k_2019.mat','sm_RamseyE_R_5k_2020.mat','sm_RamseyE_R_5k_2021.mat','sm_RamseyE_R_5k_2022.mat','sm_RamseyE_R_5k_2023.mat','sm_RamseyE_R_5k_2024.mat','sm_RamseyE_R_5k_2025.mat','sm_RamseyE_R_5k_2026.mat','sm_RamseyE_R_5k_2027.mat','sm_RamseyE_R_5k_2028.mat','sm_RamseyE_R_5k_2029.mat','sm_RamseyE_R_5k_2030.mat','sm_RamseyE_R_5k_2031.mat','sm_RamseyE_R_5k_2032.mat','sm_RamseyE_R_5k_2033.mat','sm_RamseyE_R_5k_2034.mat',}
%files={'sm_RamseyE_5k_R_1954.mat','sm_RamseyE_5k_R_1955.mat','sm_RamseyE_5k_R_1956.mat','sm_RamseyE_5k_R_1957.mat','sm_RamseyE_5k_R_1958.mat','sm_RamseyE_5k_R_1959.mat','sm_RamseyE_5k_R_1960.mat','sm_RamseyE_5k_R_1961.mat','sm_RamseyE_5k_R_1962.mat','sm_RamseyE_5k_R_1963.mat','sm_RamseyE_5k_R_1999.mat','sm_RamseyE_R_5k_2003.mat','sm_RamseyE_R_5k_2004.mat','sm_RamseyE_R_5k_2005.mat','sm_RamseyE_R_5k_2006.mat','sm_RamseyE_R_5k_2007.mat','sm_RamseyE_R_5k_2008.mat','sm_RamseyE_R_5k_2009.mat','sm_RamseyE_R_5k_2010.mat','sm_RamseyE_R_5k_2011.mat','sm_RamseyE_R_5k_2012.mat','sm_RamseyE_R_5k_2013.mat','sm_RamseyE_R_5k_2014.mat','sm_RamseyE_R_5k_2015.mat','sm_RamseyE_R_5k_2016.mat','sm_RamseyE_R_5k_2017.mat','sm_RamseyE_R_5k_2018.mat','sm_RamseyE_R_5k_2019.mat','sm_RamseyE_R_5k_2020.mat','sm_RamseyE_R_5k_2021.mat','sm_RamseyE_R_5k_2022.mat','sm_RamseyE_R_5k_2023.mat','sm_RamseyE_R_5k_2024.mat','sm_RamseyE_R_5k_2025.mat','sm_RamseyE_R_5k_2026.mat','sm_RamseyE_R_5k_2027.mat','sm_RamseyE_R_5k_2028.mat','sm_RamseyE_R_5k_2029.mat','sm_RamseyE_R_5k_2030.mat','sm_RamseyE_R_5k_2031.mat','sm_RamseyE_R_5k_2032.mat','sm_RamseyE_R_5k_2033.mat','sm_RamseyE_R_5k_2034.mat'};
%files={'sm_RamseyE_5k_R_1954.mat','sm_RamseyE_5k_R_1955.mat','sm_RamseyE_5k_R_1956.mat','sm_RamseyE_5k_R_1957.mat','sm_RamseyE_5k_R_1958.mat','sm_RamseyE_5k_R_1959.mat','sm_RamseyE_5k_R_1960.mat','sm_RamseyE_5k_R_1961.mat','sm_RamseyE_5k_R_1962.mat','sm_RamseyE_5k_R_1963.mat','sm_RamseyE_5k_R_1999.mat','sm_RamseyE_R_5k_2003.mat','sm_RamseyE_R_5k_2004.mat','sm_RamseyE_R_5k_2005.mat','sm_RamseyE_R_5k_2006.mat','sm_RamseyE_R_5k_2007.mat','sm_RamseyE_R_5k_2008.mat','sm_RamseyE_R_5k_2009.mat','sm_RamseyE_R_5k_2010.mat','sm_RamseyE_R_5k_2011.mat','sm_RamseyE_R_5k_2012.mat','sm_RamseyE_R_5k_2013.mat','sm_RamseyE_R_5k_2014.mat','sm_RamseyE_R_5k_2015.mat','sm_RamseyE_R_5k_2016.mat','sm_RamseyE_R_5k_2017.mat','sm_RamseyE_R_5k_2018.mat','sm_RamseyE_R_5k_2019.mat','sm_RamseyE_R_5k_2020.mat','sm_RamseyE_R_5k_2021.mat','sm_RamseyE_R_5k_2022.mat','sm_RamseyE_R_5k_2023.mat','sm_RamseyE_R_5k_2024.mat','sm_RamseyE_R_5k_2025.mat','sm_RamseyE_R_5k_2026.mat','sm_RamseyE_R_5k_2027.mat','sm_RamseyE_R_5k_2028.mat','sm_RamseyE_R_5k_2029.mat','sm_RamseyE_R_5k_2030.mat','sm_RamseyE_R_5k_2031.mat','sm_RamseyE_R_5k_2032.mat','sm_RamseyE_R_5k_2033.mat','sm_RamseyE_R_5k_2034.mat','sm_RamseyE_R_5k_1968.mat','sm_RamseyE_R_5k_1969.mat','sm_RamseyE_R_5k_1970.mat','sm_RamseyE_R_5k_1971.mat','sm_RamseyE_R_5k_1972.mat','sm_RamseyE_R_5k_1973.mat','sm_RamseyE_R_5k_1974.mat','sm_RamseyE_R_5k_1975.mat','sm_RamseyE_R_5k_1976.mat','sm_RamseyE_R_5k_1977.mat','sm_RamseyE_R_5k_1978.mat','sm_RamseyE_R_5k_1979.mat','sm_RamseyE_R_5k_1980.mat'}

files = uigetfile('sm*.mat', 'MultiSelect', 'on');
if ~iscell(files)
    files = {files};
end
files = sort(files);
['{' sprintf('''%s'',', files{:}),'}']
frames={};
frames{length(files)}=[];
frames{1}=1:20;

%% analyze all the data (this is slow);
clear ta2;
for i=1:length(files)
  ta2(i) = mfit_echo_fish(files{i},struct('mfitopts','','frames',frames{i}));
end
%%
taa=ta2;
for i=1:length(files)
    %taa(i).s=[];
    taa(i).mdata=[];
    taa(i).scantime=ta2(i).s.scantime;
    taa(i).mmod=[];
end
save('ana_tdep_results_2012_09_24','taa');
% Load the temp log and assign it to the data
%%
load('ana_tdep_results_2012_09_24');
%%
f=fopen('z:/qDots/fridge/matlab_fridge_log/MCLog_2011_12_19.txt');
tlog=[];
ind = 0;
while 1  tl = fgets(f);
  if ~ischar(tl);
      break;
  end
  if ind >1
  d=regexp(tl,'\s*','split');
  if ~isempty(d)
    tmp=sscanf(d{3},'%f');
    if ~isempty(tmp)
      tm=datenum(sprintf('%s %s',d{1},d{2}));
      tlog(end+1,:)=[tm ,1e3*tmp];
    end
  end
  end
  ind = ind+1;
end
fclose(f);

%% Make non-hybrid plots
% p(1)=t2
% p(2)=alpha
% p(3)=a
% p(4+5i)=w_i
% p(5+5i)=phi_i
% p(6+5i)=t2s
% p(7+5i)=y0_i
% p(8+5i)=x0_i
ta=ta2;
for i=1:length(ta)
    
    [mtd,ind] = min(abs(tlog(:,1)-ta(i).s.scantime));    
    fprintf('Temp was %g\m',mtd);
    ta(i).temp = tlog(ind,2);    
    ta(i).std=sqrt(diag(ta(i).cov))';   
    aa=ta(i).p;
    as=ta(i).std;
    ind=6:5:size(aa,2);
    %ind=[21:5:70];
    ta(i).t2s=aa(5);
    ta(i).t2s_s=as(5);
    ta(i).alpha=ta(i).p(2);
    ta(i).alphas=sqrt(ta(i).cov(2,2));
    ta(i).t2=abs(ta(i).p(1));
    ta(i).t2_s=sqrt(ta(i).cov(1,1));
end
if 0
 gind=find(([ta.dbz] < 45) & ([ta.t2] < 15) & ([ta.t2] > 1e-3) & ([ta.chisq] < 1)) ;
 find(([ta.dbz] < 45) & ([ta.t2] < 15) & ([ta.t2] > 1e-3) & ([ta.chisq] < 5))
 data=ta(gind);
else
 data=ta(find([ta.dbz] < 45 & [ta.chisq] < 5 ));
end
% plot dbz crap

opts={'MarkerSize',10,'LineWidth',2};
figure(1); clf;
plot([data.temp], [data.dbzt2], 'rx',opts{:})
%x=linspace(min([data.temp]),max([data.temp]),100);
%hold on; plot(x,1e4*x.^-2);
%hold on; plot(x,1e2*x.^-1,'r-');
xlabel('T(mK)');
ylabel('\Delta B_z T_2^* (ns)');
set(gca,'YLim',[0 120])

% plot t2s crap

opts={'MarkerSize',10,'LineStyle','None'};
figure(2); clf;
errorbar([data.temp], [data.t2s], [data.t2s_s],opts{:})
%plot([data.temp],[data.t2s],'rx')
%x=linspace(min([data.temp]),max([data.temp]),100);
%hold on; plot(x,1e4*x.^-2);
%hold on; plot(x,1e2*x.^-1,'r-');
xlabel('T(mK)');
ylabel('T_2^* extracted from echo envvelope(ns)');
set(gca,'YLim',[0 100])
%
% Temperature plot
fignum = 50;
figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
pf=@loglog;
pf(x,s*8e3*x.^-2,'g-');
hold on;
errorbar([data.temp],[data.t2],[data.t2_s],'LineStyle','None');
xlabel('T (mK)');
ylabel('T_2 (\mus)');

gind=find([data.alphas] < 5e100);
figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
errorbar([data(gind).temp],[data(gind).alpha],[data(gind).alphas],'b.','LineStyle','None');
xlabel('T (mK)');
ylabel('\beta');

figure(fignum); clf; fignum = fignum+1;
errorbar([data.temp], [data.t2s], [data.t2s_s],'b.',opts{:})
hold on;
errorbar([data.temp], [data.dbzt2], [data.dbzt2s],'r.',opts{:})
hold on;
xlabel('T(mK)');
ylabel('T_2^* extracted from echo envvelope(ns)');
set(gca,'YLim',[0 130])
%

%% Make beta-plot w/ averaging
temps=unique(sort([data.temp]));
length(data)
length(temps)
for i=1:length(temps)
    inds=find([data.temp] == temps(i));
    beta(i)=nanmean([data(inds).alpha]);
    betas(i)=sqrt(sum(([data(inds).alphas]).^2))/length(inds);
end
figure(51);
hold on;
errorbar(temps,beta,betas,'r.');
ind=find(temps < 100);
[p chisq cov]=mfitwrap(struct('x',temps(ind),'y',beta(ind),'vy',betas(ind).^2,'rng',[150 inf]), ...
      struct('fn',@(p,x) p(1)*x+p(2)),[ 0 0],'plinit plfit err')
%% Make beta-plot w/ extra averaging
tempin=unique(sort([data.temp]));
for i=1:length(tempin)  
    if tempin(i) > 100
      inds=find(abs([data.temp] - tempin(i)) < 20);
    else
      inds=find(abs([data.temp] - tempin(i)) < 1e-6);        
    end
    beta(i)=sum([data(inds).alpha]./([data(inds).alphas]))/sum(1./([data(inds).alphas]))    
    betas(i)=sqrt(length(inds))/sum(1./([data(inds).alphas]));
    temps(i)=mean([data(inds).temp]);
end
figure(51);
hold on;
errorbar(temps,beta,betas,'r.');
ind=find(temps < 100);
[p chisq cov]=mfitwrap(struct('x',temps(ind),'y',beta(ind),'vy',betas(ind).^2,'rng',[150 inf]), ...
      struct('fn',@(p,x) p(1)*x+p(2)),[ 0 0],'plinit plfit err')
  
%% Make hybrid plots
% p(1)=t2
% p(2)=t2f
% p(3)=a
% p(4+5i)=w_i
% p(5+5i)=phi_i
% p(6+5i)=t2s
% p(7+5i)=y0_i
% p(8+5i)=x0_i
ta=ta2;
for i=1:length(ta)
    
    [mtd,ind] = min(abs(tlog(:,1)-ta(i).s.scantime));    
    ta(i).temp = tlog(ind,2);    
    ta(i).std=sqrt(diag(ta(i).cov))';   
    aa=ta(i).p;
    as=ta(i).std;
    ind=6:5:size(aa,2);
    %ind=[21:5:70];
    ta(i).t2s=aa(5);
    ta(i).t2s_s=as(5);
    ta(i).alpha=ta(i).p(2);
    ta(i).alphas=sqrt(ta(i).cov(2,2));
    ta(i).t2=abs(ta(i).p(1));
    ta(i).t2_s=sqrt(ta(i).cov(1,1));
end
gind=find(([ta.dbz] < 45) & ([ta.t2] < 15) & ([ta.t2] > 1e-3) & ([ta.chisq] < 1)) ;
gind=find(([ta.dbz] < 45) & ([ta.t2] < 15) & ([ta.t2] > 1e-3) & ([ta.chisq] < 1) & ([ta.alphas] < 5))
data=ta(gind);
% plot dbz crap

opts={'MarkerSize',10,'LineWidth',2};
figure(1); clf;
plot([data.temp], [data.dbzt2], 'rx',opts{:})
%x=linspace(min([data.temp]),max([data.temp]),100);
%hold on; plot(x,1e4*x.^-2);
%hold on; plot(x,1e2*x.^-1,'r-');
xlabel('T(mK)');
ylabel('\Delta B_z T_2^* (ns)');
set(gca,'YLim',[0 120])

% plot t2s crap

opts={'MarkerSize',10,'LineStyle','None'};
figure(2); clf;
errorbar([data.temp], [data.t2s], [data.t2s_s],opts{:})
%plot([data.temp],[data.t2s],'rx')
%x=linspace(min([data.temp]),max([data.temp]),100);
%hold on; plot(x,1e4*x.^-2);
%hold on; plot(x,1e2*x.^-1,'r-');
xlabel('T(mK)');
ylabel('T_2^* extracted from echo envvelope(ns)');
set(gca,'YLim',[0 100])
%
% Temperature plot
fignum = 50;
figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
pf=@loglog;
pf(x,s*8e3*x.^-2,'g-');
hold on;
errorbar([data.temp],[data.t2],[data.t2_s],'LineStyle','None');
xlabel('T (mK)');
ylabel('T_2 (\mus)');

gind=find([data.alphas] < 5e1000);
figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
errorbar([data(gind).temp],[data(gind).alpha],[data(gind).alphas],'LineStyle','None');
xlabel('T (mK)');
ylabel('T2 (1/f)');

figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
errorbar([data(gind).temp],[data(gind).alpha],[data(gind).alphas],'r.','LineStyle','None');
hold on;
errorbar([data(gind).temp],[data(gind).t2],[data(gind).t2_s],'b.','LineStyle','None');
xlabel('T (mK)');
ylabel('T2 (1/f), white');

figure(fignum); clf; fignum = fignum+1;
x=linspace(min([data.temp]),max([data.temp]),100);
s=1
errorbar([data(gind).temp],[data(gind).alpha],[data(gind).alphas],'r.','LineStyle','None');
hold on;
errorbar([data(gind).temp],[data(gind).t2],[data(gind).t2_s],'b.','LineStyle','None');
xlabel('T (mK)');
ylabel('1/T2 (1/f), white');

%%
% fast ramsey w/ weak t dep fastramsey_2169 and friends.
% pretty echo: 2631 w/ nice envelope.

%% 
% compare hybrid to not-hybrid
nh=load('ana_tdep_results.mat');
hy=load('ana_tdep_results_hybrid.mat');
%%
for i=1:length(nh.taa)
   pnh = nh.taa(i).p;
   ph=hy.taa(i).p;
   x=linspace(min(nh.taa(i).s.tv),max(nh.taa(i).s.tv),256);   
   figure(1);
   clf;
   plot(x,exp(-abs(x/pnh(1)).^pnh(2)),'r-');
   hold on;
   plot(x,exp(-(abs(x/ph(1)) + abs(x/ph(2)).^2)),'b-');
   pause
end
