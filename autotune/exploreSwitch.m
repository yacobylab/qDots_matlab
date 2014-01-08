function exploreSwitch(config)
%function exploreSwitch(config)
% A tool to make poking around for switches easy.
% First, take a charge scan.
if ~exist('config','var')
    config=struct();
end
config=def(config,'opts','reuse');
global ess;
global tunedata;

if ~isopt(config,'reuse')
  ess.fcs=atfinechrg('x',64);
  ess.fcs=ess.fcs{1};
end

ess.coeff=fit_plane(ess.fcs);
[ess.mx,ess.my]=meshgrid(1:size(ess.fcs,2),1:size(ess.fcs,1));
ess.psub=ess.fcs-ess.mx*ess.coeff(1)-ess.my*ess.coeff(2)-ess.coeff(3);
    
figure(64);
imagesc(tunedata.chrg.scan.loops(1).rng,tunedata.chrg.scan.loops(2).rng, ess.psub);

ess.rng=[tunedata.chrg.scan.loops(1).rng*.1;tunedata.chrg.scan.loops(2).rng*.1]';

scan=tunedata.chrg.scan;
scan.loops(2).setchan={'count'};
scan.loops(2).npoints=50;
scan.loops(1).setchan={tunedata.chrg.scan.loops(1).setchan{1},tunedata.chrg.scan.loops(2).setchan{1}};
scan.loops(1).rng=[0 1];
scan.loops(1).trafofn(1).fn=@rng_trafofn;
scan.loops(1).trafofn(1).args={1};
scan.loops(1).trafofn(2).fn=@rng_trafofn;
scan.loops(1).trafofn(2).args={2};
scan.loops(1).npoints=256;
scan.loops(1).ramptime=-20e-3;
scan.loops(2).procfn.fn.fn=@sub_slope;
scan.loops(2).procfn.fn.args={1};

ess.scan=scan;

doscan

return

% Do all the work.
function doscan
  global ess;
  updatescan;
  
  while ( 1 )
    d=smrun(ess.scan);
    if any(isnan(d{1}))
        return;
    end
  end
return


% Update the picture of where we're scanning.
function updatescan
 global ess;
 global tunedata;
 
  figure(64);
  clf; hold on;
  imagesc(tunedata.chrg.scan.loops(1).rng,tunedata.chrg.scan.loops(2).rng, ess.psub);
  axis image;
  plot(ess.rng(:,1),ess.rng(:,2),'rx-');        
  ess.axis=gca;
  set(gcf,'WindowButtonDownFcn',@mousedown);
  set(gcf,'WindowButtonUpFcn',@mouseup);
return

function mousedown(f, ev)
  global ess;
  cp=get(ess.axis,'CurrentPoint');
  ess.down=cp(1,1:2);
return

function mouseup(f, ev)
  global ess;
  cp=get(ess.axis,'CurrentPoint');
  ess.up=cp(1,1:2);
  ess.rng(1,1:2)=ess.down;
  ess.rng(2,1:2)=ess.up;  
  fprintf('up\n');
  ess.rng
  updatescan
return


function v0=rng_trafofn(x,c,c0) % loop var, ignore, x=1;y=2
 global ess;
 v0=(1-x(1))*ess.rng(1,c0) + x(1)*(ess.rng(2,c0));
return

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
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

function data=sub_slope(data,opts)
  gx=diff(data);
  gx=smooth(gx,3);
  coeff=median(gx);
  data=data-(1:length(data))'*coeff;
  data=data-mean(data);
  if opts == 1
      md=min(data);
      Md=max(data);
     data=(data-md)/(Md-md);
  else
     data=data-mean(data);
  end
return


function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return


function se=cull(data)
m = median(data(:));
s = median(abs(data(:)-m));
se = find(abs(data(:)-m) < 2*s);
m = mean(data(se));
se = find(abs(data(:)-m) < 2*s);
return