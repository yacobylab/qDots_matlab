%% Use the below to reset auto-fitting.
mask_times={};
%
%% pick file

files = uigetfile('sm*.mat','MultiSelect','on');
%files={'sm_Ramsey_fast_R_5566.mat'};
%files={};
%for j=5545:5547
%for j=5545
%    files{end+1}=sprintf('sm_Ramsey_fast_R_%d',j);
%end
if ischar(files)
    files = {files};
end
mask_times={};
%% start to ana
fignum=400;
results=[];

for i=1:length(files)
  s = load(files{i});
  scan = s.scan;
  data = squeeze(nanmean(s.data{1}));
  if 0 % prescale
      for j=1:size(data,1)
          if 0
            data(j,:)=data(j,:)-min(data(j,:));
            data(j,:)=data(j,:)/max(data(j,:));
          else
            d=sort(data(j,:));
            n=size(data,2);
            dmin=d(2);
            dmax=d(end-1);
            data(j,:)=(data(j,:)-dmin)/(dmax-dmin);
          end
      end
  end
  scantime = getscantime(scan,s.data);
  if length(scan.data.pulsegroups) >1
      gd = plsinfo('xval', scan.data.pulsegroups(1).name, [], scantime);
  else
      gd = plsinfo('xval', scan.data.pulsegroups.name, [], scantime);    
  end
  
  yv  =linspace(scan.loops(1).rng(1),scan.loops(1).rng(2), scan.loops(1).npoints);
  for j = 1:length(yv)
    yv(j) = scan.loops(1).trafofn(1).fn(yv(j));
  end
  evo = gd(2,1);
  yv = evo./yv;
  yv=yv/1e-9;
  xv = gd(3,:);
  figure(1);
  clf;
  imagesc(xv,yv,data);
  
  results(i).xv = xv;
  results(i).yv = yv;
  results(i).data = data;
  results(i).file=files{i};
    
  figure(fignum); clf; hold on;
  imagesc(xv', yv,data); axis tight;
  xlabel('\epsilon (mV)');
  ylabel('T_{evo}(ns)');
  figure(fignum+1); clf; hold on;
  %plot(data(:,end), 'r');
  plot(yv,data(:,10), 'b');
 
  while(1)
      % Try fitting one epsilon first
      j= 10; % the column that you want to fit.
      pars = [];
      figure(fignum+2); clf; hold on;
      plot(data(:,j));
      if i > length(mask_times) || any(isnan(mask_times{i}))
        m = ginput(2);        
        %m=[-inf inf;inf inf];
        m(:,1)=min(max(1,m(:,1)),size(data,1));        
        mask_times{i}=round(m(1,1)):round(m(2,1));
      end
      mask=mask_times{i};

      figure(fignum+3); clf; hold on;
      params = [];
      dy=mean(range(data))/2;
      %opts='plinit plfit pause fine';
      opts='samefig woff plfit fine';
      if 0
          fitfn = @(p,x)p(1)+ exp(-(x/p(2)).^2).*(p(3)*cos(x*p(4))+p(5)*sin(p(4)*x));
          for j = 1:size(data,2)
              beta0 = fioscill(yv(mask), j*dy+data(mask,j)', 1);
              beta0(2)=range(yv(mask))/4;              
              params(j,:) = fitwrap(opts, yv(mask),j*dy+data(mask,j)', beta0, fitfn);
          end
           results(i).freq=1e3*params(:,4)/(2*pi);
      elseif 0  % fit chirp
          fitfn = @(p,x)p(1)+ exp(-(x/p(2)).^2).*(p(3)*cos(x.*(p(4)+p(6)*exp(-x/p(7))))+p(5)*sin(x.*(p(4)+p(6)*exp(-x/p(7)))));
          for j = 1:size(data,2)
              beta0 = fioscill(yv(mask), j*dy+data(mask,j)', 1);
              beta0(6)=beta0(4)/3;              
              beta0(2)=range(yv(mask));
              beta0(7)=beta0(2)
              params(j,:) = fitwrap(opts, yv(mask),j*dy+data(mask,j)', beta0, fitfn);
              params(j,:)              
              %input('Hi','s');
          end
           results(i).freq=1e3*params(:,4)/(2*pi);
      elseif 0 % Analyze manually.
        figure(401); clf;
        for j=1:size(data,2)
          plot(yv(mask),data(mask,j))
          d=diff(ginput(2))
          results(i).freq(j)=1/d(2);
        end            
      else % Analyze a data set by overlapping psd.
        fs=1/(yv(2)-yv(1));
        h=spectrum.welch('Hamming',64);
        %h=spectrum.music(2);
        figure(401);      clf;        
        figure(402);      clf;        
        img=[];
        %fitfunc=@(p,x) p(1)+p(2)./(((x-p(3))/p(4)).^2+1)        
        %fitfunc=@(p,x) p(1)+p(2)*exp(-(((x-p(3))/p(4)).^2)        )
        fitfunc=@(p,x) p(1)+p(2)*exp(-(((x-p(3))/p(4)).^2))+p(5)*exp(-(x/p(6).^2));
        %pts=find(xv > 0 & xv < 1.5);
        pts=1:size(data,2)        
        for j=pts
            hpsd{j}=psd(h,data(mask,j)-mean(data(mask,j)),'Fs',fs);
            %hpsd{j}=pseudospectrum(h,data(mask,j)-mean(data(mask,j)),'Fs',fs);
            figure(401);
            plot(hpsd{j});
            hold on;            
            img(:,j)=hpsd{j}.Data; 
            skip=5;
            y=hpsd{j}.Data(skip:end) + j*1e-3;
            x=hpsd{j}.Frequencies(skip:end);
            if 0
              figure(402);    
              clf;
              plot(x,log(y));
              g=ginput(1);
            else
              [m mi] = max(y);
              g(1)=x(mi);
            end
            params(j,:)=fitwrap('plot plinit plfit',x',y', [min(y) (max(y)-min(y)) g(1) 1 (y(2)-min(y)) 2],fitfunc);
            hold on;                          
            if params(j,4) > range(x)/2 || params(j,3) < 0 || params(j,3) > max(x)
                params(j,3:4)=nan;
            end
        end
        for j=pts(end)+1:size(data,2)
            params(j,:)=nan;
        end
        figure(403);
        clf;
        results(i).freq=1e9*params(:,3);
        imagesc(xv,hpsd{1}.Frequencies,log(img));
        colormap(bone);
        hold on;
        errorbar(xv,params(:,3), params(:,4),'r-');
        axis tight;
      end
      
      figure(fignum+4); clf; hold on;
      plot(xv, 1e3*abs(params(:,4))/(2*pi),'b-x');
      set(gca,'YLim',[0 40000]);
      xlabel('\epsilon(mV)'); ylabel('J(MHz)');
      if strcmp(input('Accept (y/N)','s'),'y')          
          break;
      else
          mask_times{i}=nan;
      end
  end
end 
%% start to ana (this is for new scans)
fignum=400;
results=[];

for i=1:length(files)
  s = load(files{i});
  scan = s.scan;
  data = squeeze(nanmean(s.data{1}));
  if 0 % prescale
      for j=1:size(data,1)
          if 0
            data(j,:)=data(j,:)-min(data(j,:));
            data(j,:)=data(j,:)/max(data(j,:));
          else
            d=sort(data(j,:));
            n=size(data,2);
            dmin=d(2);
            dmax=d(end-1);
            data(j,:)=(data(j,:)-dmin)/(dmax-dmin);
          end
      end
  end
  scantime = getscantime(scan,s.data);
  if length(scan.data.pulsegroups) >1
      gd = plsinfo('xval', scan.data.pulsegroups(2).name, [], scantime);
  else
      gd = plsinfo('xval', scan.data.pulsegroups.name, [], scantime);    
  end
  
%   yv  =linspace(scan.data.frange(1),scan.data.frange(2), scan.loops(1).npoints);
%   for j = 1:length(yv)
%     yv(j) = scan.loops(1).trafofn(1).fn(yv(j));
%   end

  yv =scan.data.freqs;
  evo = gd(2,1);
  yv = evo./yv;
  yv=yv/1e-9;
  xv = gd(3,:);
  figure(1);
  clf;
  imagesc(xv,yv,data);
  
  results(i).xv = xv;
  results(i).yv = yv;
  results(i).data = data;
  results(i).file=files{i};
    
  figure(fignum); clf; hold on;
  imagesc(xv', yv,data); axis tight;
  xlabel('\epsilon (mV)');
  ylabel('T_{evo}(ns)');
  figure(fignum+1); clf; hold on;
  %plot(data(:,end), 'r');
  plot(yv,data(:,10), 'b');
 
  while(1)
      % Try fitting one epsilon first
      j= 10; % the column that you want to fit.
      pars = [];
      figure(fignum+2); clf; hold on;
      plot(data(:,j));
      if i > length(mask_times) || any(isnan(mask_times{i}))
        m = ginput(2);        
        %m=[-inf inf;inf inf];
        m(:,1)=min(max(1,m(:,1)),size(data,1));        
        mask_times{i}=round(m(1,1)):round(m(2,1));
      end
      mask=mask_times{i};

      figure(fignum+3); clf; hold on;
      params = [];
      dy=mean(range(data))/2;
      %opts='plinit plfit pause fine';
      opts='samefig woff plfit fine';
      if 0
          fitfn = @(p,x)p(1)+ exp(-(x/p(2)).^2).*(p(3)*cos(x*p(4))+p(5)*sin(p(4)*x));
          for j = 1:size(data,2)
              beta0 = fioscill(yv(mask), j*dy+data(mask,j)', 1);
              beta0(2)=range(yv(mask))/4;              
              params(j,:) = fitwrap(opts, yv(mask),j*dy+data(mask,j)', beta0, fitfn);
          end
           results(i).freq=1e3*params(:,4)/(2*pi);
      elseif 0  % fit chirp
          fitfn = @(p,x)p(1)+ exp(-(x/p(2)).^2).*(p(3)*cos(x.*(p(4)+p(6)*exp(-x/p(7))))+p(5)*sin(x.*(p(4)+p(6)*exp(-x/p(7)))));
          for j = 1:size(data,2)
              beta0 = fioscill(yv(mask), j*dy+data(mask,j)', 1);
              beta0(6)=beta0(4)/3;              
              beta0(2)=range(yv(mask));
              beta0(7)=beta0(2)
              params(j,:) = fitwrap(opts, yv(mask),j*dy+data(mask,j)', beta0, fitfn);
              params(j,:)              
              %input('Hi','s');
          end
           results(i).freq=1e3*params(:,4)/(2*pi);
      elseif 0 % Analyze manually.
        figure(401); clf;
        for j=1:size(data,2)
          plot(yv(mask),data(mask,j))
          d=diff(ginput(2))
          results(i).freq(j)=1/d(2);
        end            
      else % Analyze a data set by overlapping psd.
        fs=1/(yv(2)-yv(1));
        h=spectrum.welch('Hamming',64);
        %h=spectrum.music(2);
        figure(401);      clf;        
        figure(402);      clf;        
        img=[];
        %fitfunc=@(p,x) p(1)+p(2)./(((x-p(3))/p(4)).^2+1)        
        %fitfunc=@(p,x) p(1)+p(2)*exp(-(((x-p(3))/p(4)).^2)        )
        fitfunc=@(p,x) p(1)+p(2)*exp(-(((x-p(3))/p(4)).^2))+p(5)*exp(-(x/p(6).^2));
        %pts=find(xv > 0 & xv < 1.5);
        pts=1:size(data,2);        
        for j=pts
            hpsd{j}=psd(h,data(mask,j)-mean(data(mask,j)),'Fs',fs);
            %hpsd{j}=pseudospectrum(h,data(mask,j)-mean(data(mask,j)),'Fs',fs);
            figure(401);
            plot(hpsd{j});
            hold on;            
            img(:,j)=hpsd{j}.Data; 
            skip=5;
            y=hpsd{j}.Data(skip:end) + j*1e-3;
            x=hpsd{j}.Frequencies(skip:end);
            if 0
              figure(402);    
              clf;
              plot(x,log(y));
              g=ginput(1);
            else
              [m mi] = max(y);
              g(1)=x(mi);
            end
            params(j,:)=fitwrap('plot plinit plfit',x',y', [min(y) (max(y)-min(y)) g(1) 1 (y(2)-min(y)) 2],fitfunc);
            hold on;                          
            if params(j,4) > range(x)/2 || params(j,3) < 0 || params(j,3) > max(x)
                params(j,3:4)=nan;
            end
        end
        for j=pts(end)+1:size(data,2)
            params(j,:)=nan;
        end
        figure(403);
        clf;
        results(i).freq=1e9*params(:,3);
        imagesc(xv,hpsd{1}.Frequencies,log(img));
        colormap(bone);
        hold on;
        errorbar(xv,params(:,3), params(:,4),'r-');
        axis tight;
      end
      
      figure(fignum+4); clf; hold on;
      plot(xv, 1e3*abs(params(:,4))/(2*pi),'b-x');
      set(gca,'YLim',[0 40000]);
      xlabel('\epsilon(mV)'); ylabel('J(MHz)');
      if strcmp(input('Accept (y/N)','s'),'y')          
          break;
      else
          mask_times{i}=nan;
      end
  end
end    
%%
figure(1);
clf;
colors='rgbcmyk';
scales=ones(1,length(results));
scales(3)=3;
scales(4)=3;
for    i = 1:length(results)
      plot(results(i).xv,scales(i)*results(i).freq,[colors(1+mod(i-1,end)) 'x-']);
  hold on;
end
x=get(gca,'XLim');
x=linspace(x(1),x(2),2000);
plot(x,(j(x) < 10000) .* j(x),'r-');
l=legend(results.file);
set(l,'interpreter','none');
%%
s = load(uigetfile('sm*.mat'));

scan = s.scan;
data = squeeze(mean(s.data{1}));
%%
scantime = getscantime(scan,s.data);
gd = plsinfo('xval', scan.data.pulsegroups.name, [], scantime);
yv  =linspace(scan.loops(1).rng(1),scan.loops(1).rng(2), scan.loops(1).npoints);
for j = 1:length(yv)
    yv(j) = scan.loops(1).trafofn(1).fn(yv(j));
end
evo = gd(2,1);
yv = evo./yv;
xv = gd(3,:);

fignum = 432;
%%
figure(fignum); clf; hold on;
imagesc(xv', yv,data); axis tight;
xlabel('\epsilon (mV)');
ylabel('T_{evo}(ns)');
figure(fignum+1); clf; hold on;
%plot(data(:,end), 'r');
plot(yv,data(:,10), 'b');


%% Try fitting one epsilon first
j= 10; % the column that you want to fit. 

pars = [];
figure(fignum+2); clf; hold on; 
plot(data(:,j));
m = ginput(2); 
if m(1,1)<0
    m(1,1) = 0;
end
if m(2,1) > size(data,1)
    m(2,1)=size(data,1);
end
mask = (round(m(1,1)):round(m(2,1))); 
fitfn = @(p,x)p(1)+ exp(-x/p(2)).*(p(3)*cos(x*p(4))+p(5)*sin(p(4)*x));
 %:size(data,1)
 beta0 = fioscill(yv(mask), data(mask,j)', 1);
 %params = fitwrap('plinit plfit', yv(end/2:end)',data(end/2:end,j), beta0, fitfn, [1 1 1 1]);
 pars = fitwrap('woff plinit plfit', yv(mask),data(mask,j)', beta0, fitfn);
 
 %% now fit all of them 
% params = [];
% figure(fignum+2); clf; hold on; 
% plot(data(:,1));
% m = ginput(1); mask = (round(m(1)):length(yv)); 
fitfn = @(p,x)p(1)+ exp(-x/p(2)).*(p(3)*cos(x*p(4))+p(5)*sin(p(4)*x));
params = [];
for j = 1:size(data,2)
beta0 = fioscill(yv(mask), data(mask,j)', 1);
params(end+1,:) = fitwrap('woff plinit plfit', yv(mask),data(mask,j)', beta0, fitfn);
end

figure(fignum+3); clf; hold on;
plot(xv, 1e-6*params(:,4)/(2*pi));
xlabel('\epsilon(mV)'); ylabel('J(MHz');
%% keep track of them
dummy = [1e-6*params(:,4)/(2*pi), xv'];
masterP = [masterP;dummy];
%%
figure(fignum+4); clf; hold on;
plot(masterP(:,7), masterP(:,4)/(2*pi), '.b');
xlabel('\epsilon(mV)'); ylabel('J(MHz');

%% now get and fit some regular ramsey data;
[figs pars]= ana_echo('',struct('opts', 'ramsey frq noppt', 't1', 1/50));
pars(:,:) = pars(:, [2,1]);
masterP = [masterP; pars'];
%% plot everything and fit it
figure(1); clf; hold on;
sets = [{1:123}, {124:size(masterP,1)}];
for i= 1:length(sets)
plot(masterP(sets{i},2), masterP(sets{i},1), [c(i) '.']);
end

beta = [mean(masterP(:,1)), .5*range(real(masterP(:,1))),.5*range(masterP(:,2))]; 
ffn = @(p,x) p(1)+p(2)*exp(-x/p(3));
params = fitwrap('plinit plfit', masterP(1:123,2)', masterP(1:123,1)',beta, ffn);

figure(1);
plot(masterP(:,2), ffn(params, masterP(:,2)));
xlabel('\epsilon(mV)'); ylabel('J(Hz');
%%
c = 'rgbcmyk';
figure(fignum+4); clf; hold on;
scales=[1 1 4 4 4 4 4 4 4 4];
for i = 1:41:size(masterP,1)
plot(masterP(i:i+40,7), scales(1+(i-1)/41)*masterP(i:i+40,4)/(2*pi), [c(1+mod(i,length(c))) '.']);
hold on;
end
xlabel('\epsilon(mV)'); ylabel('J(Hz');
legend('1','2','3','4','5','6','7','8','9');

%%
figure(1); clf; hold on;
sets = [{1:7}, {8:14}, {15:82}];

for i =1:2%length(sets)
    plot(p(1,sets{i}), p(2, sets{i}), [c(i) '.']);
end


%% extract T2 from lots of RamseyE scan;
files = [];
T2data = [];
for j=5554:5559
    files{end+1}=sprintf('sm_RamseyE_R_%d',j);
end

for i = 1:length(files)
    [f p] = ana_echo(files{i},struct('opts', 'noppt amp per', 't1', 1/50, 'grng', [.008 inf]));
    s = load(files{i});
    scantime = getscantime(s.scan,s.data);
    %gd = plsinfo('xval', s.scan.data.pulsegroups(2).name, [], scantime);
    gd2 = plsinfo('params', s.scan.data.pulsegroups(2).name, [], scantime);
    eps = gd2(1);
    T2data(i).eps = eps;
    T2data(i).T2 = p.T2;
    T2data(i).file = files{i};    
end
fnum = 684;
figure(fnum);
p = plot([T2data(:).eps], [T2data(:).T2], 'b.');
set(p,'MarkerSize', 10);

