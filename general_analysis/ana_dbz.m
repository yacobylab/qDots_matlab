function val = ana_dbz(file, config)

if ~exist('file','var') || isempty(file)
    file = uigetfile('sm*.mat');
end

if ~exist('config','var')
    config=struct();
end

fprintf('%s \n', file);
%load file
s=load(file);

config = def(config,'opts','fid');   % Random boolean options
config = def(config,'t1Time',{'before'}); % defaults to try to get t1 from before 
%populate all of the t1times
while length(config.t1Time)<length(s.scan.loops(1).getchan)-1
    config.t1Time = [config.t1Time config.t1Time{end}];
end
scantime=getscantime(s.scan,s.data);
data = {};
figure(44); clf;

for j = 1:length(s.scan.loops(1).getchan)-1

  if strcmp(s.scan.loops(1).getchan{j}(end), '1')
      side = 'left';
  elseif strcmp(s.scan.loops(1).getchan{j}(end), '2')
      side = 'right';
  else
      error('could not determine proper datachans');
      return;
  end
  
[t1t t1] = att1(side,scantime,config.t1Time{j});
d=anaHistScale(s.scan,{s.data{j}},t1);
data{j} = squeeze(mean(d{1},1));

subplot(length(s.scan.loops(1).getchan)-1,1,j); hold on;
plot(data{j});

end


end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
return;
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end