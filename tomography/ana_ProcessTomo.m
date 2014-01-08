function out = ana_ProcessTomo(filename, opts)
if ~exist('filename','var') || isempty('filename')
    filename=uigetfile('sm_processTomo*.mat');
end
if ~exist('opts','var')
    opts=struct();
end

opts=def(opts,'tomocal','');
if isempty(opts.tomocal)
    opts.tomocal=uigetfile('sm_TomoCal_rot*_mat.mat');
end
if opts.tomocal ~= 0
  tomocal=load(opts.tomocal)
else
  tomocal.data_corrector=@(x) x;
end

d=ana_avg(filename,'noplot');
figure(253);
[r1 r2 r3] = size(d.data{1});
m=reshape(permute(d.data{1},[1 3 2]), [r1,r2*r3]);
imagesc(m);
figure(254);
imagesc(squeeze(nanmean(d.data{1},1)));
colorbar
out.d=d;
axismap={'UD','Y','ST'};
for i = 1:length(d.scan.data.pulsegroups)
    for j=1:length(axismap)
       if regexp(d.scan.data.pulsegroups(i).name,axismap{j})
           fprintf('Group %d(%s) measures %s goes to %d\n', ...
              i,d.scan.data.pulsegroups(i).name,axismap{j},j);
           data(j,:) = -(squeeze(nanmean(d.data{1}(:,i,:),1))*2-1);
       end
    end
end
out.data=data;
initial=data(:,1:2:end);
final=data(:,2:2:end);
initial = tomocal.data_corrector(initial')';
final = tomocal.data_corrector(final')';
plotStates(initial,final);
return

function plotStates(initial, final)

% dbloch style
%
c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','sphere','color',[c .05])];
cmds=[cmds struct('type','label','val',[0 0 1.15 ],'label','|S>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|T>')];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|\downarrow\uparrow>')];
opts=struct('opts',' resize raise','dt',.1);
base=cmds;
opengl software;
figure(255);
clf; 
dbloch(cmds,opts);
colors='rgbcmkg';
for i=1:size(initial,2)
   hold on;
   plot3(initial(1,i),initial(2,i),initial(3,i),[colors(i) 'x']);
   plot3(final(1,i),final(2,i),final(3,i),[colors(i) 'o']);
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