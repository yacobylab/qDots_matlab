function output = atplschk2(pg,config)
% function atplschk(pg,config)
% plots pulse on top of relevant charge diagrams.
% config is a config struct with fields:
%   opts can include flat, diff, noleft, noright, markers
%   pulses is a list of pulses to plot.  defaults to 1.
%   run is a tunedata run.  Defaults to las.
%   awg is which awg to render.  Defaults to 1
%   offset is an offset to apply to the pulse (ie to check zoom pulse) defaults to zero.  1x4 vector, xl,yl,xr,yr
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

% eventually want to adapts to accept pulse number, pulsegroup, or pinf
% struct,  right now it only accepts pulsegroup. 

global tunedata; global plsdata;

if ~exist('config','var')
    config=1;
end
if isreal(config)
    config=struct('pulses',config);
end

if iscell(config)
    config=struct(config{:});
end

config=def(config,'run',[length(tunedata.sets(1).runs) length(tunedata.sets(2).runs)]);
config=def(config,'offset',[0 0 ; 0 0]);
config=def(config,'pulses',1);
config=def(config,'opts','flat');
config=def(config,'awg',1);

if isopt(config,'measp')
   config.offset = -1e3*[tunedata.sets(2).measp(1:2,1)'; tunedata.sets(1).measp(1:2,1)'];
end


if(~iscell(pg))
    pg={pg};
end

runname = {['_', tunedata.sets(1).name], ['_', tunedata.sets(2).name]}; %name of current set
%if ~isempty(runname)
%   runname{:} = ['_', runname{:}];
%end

figure(55); %clf;
set(figure(55),'Name','Pulse Check');
% nsets=2;
% for j = 1:nsets
%     d=[];
%     for sj=(1:config.run(j))-1;
%       file = sprintf('%s/sm_chrg%s_%03i.mat', tunedata.dir, runname{j}, config.run(j)-sj);
%       if exist(file,'file')          
%         d = load(file);
%         break;
%       end
%     end
%     if isempty(d)
%         continue;
%     end
%     subplot(2,2,1+nsets-j);
%     rng(1,:)=d.scan.loops(1).rng;
%     rng(2,:)=d.scan.loops(2).rng;
%     
%     
%     %file = sprintf('%s/sm_chrg%s_%03i', tunedata.dir, runname, run);
%     %d=load(file);
%     %subplot(221);
%     %rng(1,:)=d.scan.loops(1).rng;
%     %rng(2,:)=d.scan.loops(2).rng;
%     
%     cdata=d.data{1};
%     if isopt(config,'diff')
%         cdata=diff(cdata);
%     end
%     if ~isopt(config,'noflat')
%         coeff=fit_plane(cdata);
%         [mx,my]=meshgrid(1:size(cdata,2),1:size(cdata,1));
%         cdata=cdata-mx*coeff(1)-my*coeff(2)-coeff(3);
%     end
%     
%     imagesc(rng(1,:)*1e3+config.offset(1+nsets-j,1),rng(2,:)*1e3+config.offset(1+nsets-j,2),cdata);
%     axis image;
%     set(gca,'ydir','normal');
%     hold on;
%     
% end

% Flatten the pulsegroups by substituting in for sequence combined groups.
changed = 1;
pg2=pg;
for k=1:length(pg2)
  pg{k}.name=pg2{k};
  pg{k}.pulses=config.pulses;
end
while changed
   for k=1:length(pg)
   changed = 0;
    if ~ischar(pg{k}.name)
        global awgdata;
        fprintf('Pulsegroup %d is %s\n',pg{k}.name,awgdata.pulsegroups(pg{k}).name);
        pg{k}.name=awgdata.pulsegroups(pg{k}).name;
        changed = 1;
    end    
   
    dd=plsinfo('gd',pg{k}.name);
    if ~isempty(strfind(dd.ctrl,'seq'))
        lg=length(dd.pulses.groups);
        pgsrc=pg{k};
        pg{k+lg-1:end+lg-1}=pg{k:end};        
        for l=1:lg
          fprintf('Hybridizing in %s\n',dd.pulses.groups{l});
          pg{k+l-1}.name = dd.pulses.groups{l};
          pg{k+l-1}.pulses=dd.pulseind(l,pgsrc.pulses);
        end
        changed = 1;
        continue
    end
   end
end

% Load the pulses
c1=[1 0 0]; % Start and end color for time traces.
c2=[0 1 0];
for k=1:length(pg)
    dd=plsmakegrp(pg{k}.name,'',pg{k}.pulses);
        
    for i=1:length(dd.pulses)
        if(isfield(dd,'pulseind'))
          d=plstowf(dd.pulses(dd.pulseind(i)));
        else
          d=plstowf(dd.pulses(i));
        end
        outchans=cell(4,1);
        outmark=cell(4,1);
        if isopt(config,'noleft')
            cmin=3;
        else
            cmin=1;
        end
        if isopt(config,'noright')
            cmax=2;
        else
            cmax=inf;
        end
        awg=config.awg;
        for l=max(1,cmin):min(length(dd.chan),cmax)
            outchans{dd.chan(l)}=d.data(awg).wf(l,:);
            outmark{dd.chan(l)}=d.data(awg).marker(l,:);
        end        
%         if(~isempty(outchans{2}) && ~isempty(outchans{1}) && ~isopt(config,'noleft'));
%             subplot(221);
%             pd=plot(outchans{2},outchans{1});
%             if isopt(config,'markers')
%                 plot(outchans{2},outchans{1},'kx');
%             end
%             
%             colorplot(outchans{2},outchans{1},c1,c2);
%         end
%         if(~isempty(outchans{3}) && ~isempty(outchans{4}) && ~isopt(config,'noright'))
%             subplot(222);
%             if isopt(config,'markers')
%                 plot(outchans{3},outchans{4},'kx');
%             end
%             colorplot(outchans{3},outchans{4},c1,c2);
%         end
        
        colors={[1 0 0],[0.5 0 0], [0 1 0], [0 0.5 0]};        
        hold on;
        for l=1:length(outchans)
          p=plot(outchans{l});
          hold on;
          p2 = plot(outmark{l});
          set(p,'Color',colors{l});
          set(p2,'Color',colors{l});
        end
%        
    end
end
output=[outchans,outmark];
return;
function pd=colorplot(x,y,c1,c2)
  ns=10;
  l=round(linspace(1,length(x),ns));
  c=linspace(0,1,ns-1);
  for i=2:ns
      pd(i)=plot(x(l(i-1):l(i)),y(l(i-1):l(i)));
      set(pd(i),'Color',c1*(1-c(i-1))+c2*c(i-1));
  end
return
function coeff=fit_plane(data)
[gx,gy] = gradient(data);
sm=2;
for l=1:size(gx,1)
    gx(l,:)=smooth(gx(l,:),3);
end
for l=1:size(gy,2)
    gy(:,l)=smooth(gy(:,l),3);
end
coeff(1)=median(gx(cull(gx)));
coeff(2)=median(gy(cull(gy)));
coeff(3)=mean(mean(data));
return

function se=cull(data)
m = median(data(:));
s = median(abs(data(:)-m));
se = find(abs(data(:)-m) < 2*s);
m = mean(data(se));
se = find(abs(data(:)-m) < 2*s);
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
