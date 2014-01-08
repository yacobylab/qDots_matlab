%% get the files

files = uigetfile('sm*.mat', 'MultiSelect', 'on');
if ~iscell(files)
    files = {files};
end
files = sort(files);
['{' sprintf('''%s'',', files{:}),'}'];
%% analyze all the data (this is slow);

results = [];
alldata=[];
grng={[1 Inf]};%{[.5 1.5]};
for j = [1:5 7]%[1:5 7:length(files)];%[length(files):-4:1 length(files)-1:-4:1]
%  try
    [f p sd]=ana_echo(files{j},struct('opts', 'echo noppt power'));  
    %[f p sd]=ana_echo(files{j},struct('opts', 'echo noppt', 't1', 1e-3));  
    %[f p sd]=ana_echo(files{j},struct('opts', 'echo noppt logfit powerlaw fitdecay rmoutlier fitoffset', 't1', 1e-3, 'grng', grng{min(j,end)}));

    results{end+1} = p;
    alldata{end+1}=sd;
    %pause
%  catch err
%      results{end+1} = [];
%      fprintf('Fit error: %s\n',err.message);
%      keyboard
%  end
end
close all;


%% Load the temp log and assign it to the data

%f=fopen('z:/qDots/fridge/matlab_fridge_log/MCLog_2011_12_19.txt');
f=fopen('z:/qDots/fridge/matlab_fridge_log/MCLog_2012_10_03.txt');
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
data=[];
%is = ~
%
for i=1:length(results)
    
    data(i).scantime = results{i}.scantime;
    [mtd,ind] = min(abs(tlog(:,1)-data(i).scantime));    
    data(i).temp = tlog(ind,2);    
    data(i).q = results{i}.Q;
    data(i).exponent = results{i}.afp(end);
    data(i).dbzt2=results{i}.dbzt2; 
    data(i).fnr = sscanf(files{i}(end-7:end-4),'%d');  
    data(i).Jbar = results{i}.Jbar;
    data(i).t2 = results{i}.T2;
    data(i).params = results{i}.params;
    data(i).t2s = mean(results{i}.params(2,6));
    data(i).scanpars=alldata{i}.scan.data;
    if isfield(data(i).scanpars,'RF1')
        data(i).RF1=data(i).scanpars.RF1;
    end
end
%%
figind = 50;
opts={'MarkerSize',10,'LineWidth',2};
figure(figind); clf; hold on;
plot([data.temp], [data.t2], 'rx',opts{:})
%x=linspace(min([data.temp]),max([data.temp]),100);
fitfn = @(p,x) p(1)*x.^p(2)+p(3); beta0 = [data(1).t2*data(1).temp^2, -2, 0];

if 1 % allow fitting of an offset
    pars = fitwrap('plinit plfit',[data.temp],[data.t2],beta0,fitfn,[1 1 1]);
else % no offset
    pars = fitwrap('plinit plfit',[data.temp],[data.t2],beta0,fitfn,[1 1 0]);
end

figure(figind); figind = figind+1;
plot([data.temp],fitfn(pars,[data.temp]));
xlabel('T(mK)');
ylabel('T2 (\mus)');
title(sprintf('T_2 ~(Temperature)^{%.2f}',pars(2)));

figure(figind); clf; figind = figind+1;
loglog([data.temp],[data.t2],'rx',opts{:}); hold on;
loglog([data.temp],fitfn(pars,[data.temp]));
xlabel('T(mK)');
ylabel('T2 (\mus)');
title(sprintf('T_2 ~(Temperature)^{%.2f}',pars(2)));


%% beta plot
figure(figind); figind = figind+1; hold on;
plot([data.temp],[data.exponent],'rx',opts{:});
xlabel('T(mK)');
ylabel('decay exponent');

