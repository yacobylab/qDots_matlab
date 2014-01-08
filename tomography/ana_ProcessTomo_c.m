function out = ana_ProcessTomo_c(filename, opts)
% function out = ana_ProcessTomo_c(filename, opts)
% Returns out.
%   data -- uncorrected data
%   datav -- uncorrected variances
%   cdata -- corrected data
%   cdatav -- corrected variances
% opts has fields
%   tomocal: 0 for no cal, blank for ui get, filename for tomocal file.

if ~exist('filename','var') || isempty('filename')
    filename=uigetfile('sm_processTomo_rot*.mat;sm_attr_processTomo_rot*.mat;*proctomo*.mat;*ProcessTomo*.mat');
end
if ~exist('opts','var')
    opts=struct();
end
if ischar(opts)
    opts=struct('opts',opts);
end

opts=def(opts,'tomocal','');
opts=def(opts,'opts','plot');
if isempty(opts.tomocal)
    opts.tomocal=uigetfile('sm_TomoCal_rot*_mat.mat;sm_attr_TomoCal_rot*_mat.mat');
end
if opts.tomocal ~= 0
  tomocal=load(opts.tomocal);
else
  tomocal.data_corrector=@(x) x;
end
fprintf('ana_ProcessTomo_c(''%s'',struct(''tomocal'',''%s''))\n', ...
        filename, opts.tomocal);
d=ana_avg(filename,'noplot noppt');
if isopt(opts,'plot')
  figure(253);
  [r1 r2 r3] = size(d.data{1});
  m=reshape(permute(d.data{1},[1 3 2]), [r1,r2*r3]);
  imagesc(m); 
  figure(254);
  imagesc(squeeze(nanmean(d.data{1},1)));
  colorbar
end
out.d=d;
out.tv=[];
axismap={'UD','Y','ST'};
k=0;
for i = 1:length(d.scan.data.pulsegroups)
    if ~isempty(regexp(d.scan.data.pulsegroups(i).name,'^dBz'))
        continue;
    end
    for j=1:length(axismap)
       if regexp(d.scan.data.pulsegroups(i).name,axismap{j})
           if ~isopt(opts,'quiet')
             %fprintf('Group %d(%s) measures %s goes to %d\n', ...
             % i,d.scan.data.pulsegroups(i).name,axismap{j},j);
           end
           data{floor(k/3)+1}(j,:) = -(squeeze(nanmean(d.data{1}(:,i,:),1))*2-1);
           datav{floor(k/3)+1}(j,:) = (2*nanstd(d.data{1}(:,i,:),[],1)/sqrt(size(d.data{1},1))).^2;
       end
    end
    if length(d.tv) >= i
      out.tv(floor(k/3)+1) = d.tv(i);
    else
      out.tv(i) = nan;
    end
    k=k+1;
end
out.data=data;
out.datav=datav;
for i=1:length(data)
    data{i}=tomocal.data_corrector(data{i}')';
    for j=1:size(datav{i},2)
        m=diag(datav{i}(:,j));
        m=(tomocal.pm')*m*tomocal.pm;
        cdatav{i}(:,j)=diag(m);
    end
end
out.cdata=data;
out.cdatav=cdatav;

if isopt(opts,'plot')
  ana_plotTomoStates(data,out.tv);
end
return


% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;