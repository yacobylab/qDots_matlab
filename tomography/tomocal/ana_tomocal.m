function [p p_cov rd rdp data_corrector]=ana_tomocal(rtfile, config)
% function ana_tomocal(rtfile, config)
% Try to unscrambles axes based on experimental data.
% sm_TomoCal_rot_R_5281.mat is nice 3-J axis
% Axis mapping: ST=X, UD/DU=Z, Y=y, see ana_tomocal_2 for ST=Z
% Load the trivial tomography
%  config.opts -- options that control fitting.
%     ref -- replace automated singlet/triplet level fitting with a reference group
%            if available.
if ~exist('config','var') || isempty(config)
    config=struct();
end
config=def(config,'opts','ref');
config=def(config,'initial',[]);
config
if ~exist('rtfile','var') || isempty(rtfile)
  [rtfile dir]=uigetfile('sm_TomoCal_rot*.mat','Rotating Set');
  rtfile = [dir rtfile];    
end
rt=load(rtfile);
scantime=getscantime(rt.scan,rt.data);

% Make some color plots for switch checks
figure(1);
clf;
d=permute(rt.data{1},[1 3 2]);
imagesc(reshape(d,size(d,1),size(d,2)*size(d,3)));

% Scaling for non-raw calibration data
if ~isempty(strfind(rtfile,'_L_'))
    side='left';
else
    side='right';
end
[t1t t1] = att1(side,scantime,'ask' );
%[t1t t1] = att1(side,scantime,'after');
%t1 = 4.5/55;  % t1 on left.

% Exclude dbz from histograms.
dbzgrps = find(cellfun(@(p) ~isempty(p),regexp({rt.scan.data.pulsegroups.name},'^([dD][bB][zZ])')));
refgrps = find(cellfun(@(p) ~isempty(p),regexp({rt.scan.data.pulsegroups.name},'^(ref)')));
rotgrps = setdiff(1:length(rt.scan.data.pulsegroups),union(dbzgrps,refgrps));
[rtd rtff mv fp]=anaHistScale(rt.scan,rt.data,t1,rotgrps);
rtd=rtd{1}*2-1;
rtff=rtff{1};

% Definitions to help us keep axes straight. 

axismap={'ST','Y','UD'};
axisorder=[1 3 2];
%
% Try to fit the data, all basis vectors, polar coords
clear mfdata;
clear mfmodel;
clear pts;
clear rs;
clear rd;
clear srd;
minus=0; % always set to zero.
covrd={};
skip=0;
dsets=1:floor(length({rt.scan.data.pulsegroups.name})/3)-skip;
%dsets([3])=[];
%dsets(end-1:end)=[];

fprintf('%d rotation groups detected\n',floor(length({rt.scan.data.pulsegroups.name})/3));
for i=1:length(dsets)
  rs{i}=find(~cellfun(@(x)isempty(x),strfind({rt.scan.data.pulsegroups.name},sprintf('rot_%d',dsets(i)+skip))));
  rd{i}=squeeze(nanmean(rtd(:,rs{i}(axisorder),:),1))';  
  if minus
    rd{i}(:,1)=-rd{i}(:,1);
    rd{i}(:,3)=-rd{i}(:,3);
  else  
    rd{i}=-rd{i};
  end
  srd{i}=squeeze(nanstd(rtd(:,rs{i}(axisorder),:),1,1))'/sqrt(size(rtd,1));
  for j=1:size(rtd,3)
    tmp=squeeze(rtd(:,rs{i}(axisorder),j));
    covrd{i,j}=cov(tmp)/size(rtd,1);
  end
  if 0
    pts{i}=1:size(rd{i},1);
  else
    pts{i}=1:32;
  end
  mfdata(i).x=pts{1};

  xvals=plsinfo('xval',rt.scan.data.pulsegroups(rs{i}(1)).name);
  if xvals(end/2+1) == xvals(1)  % This is folded data...
    fprintf('Group %d appears to be folded\n',i);
    mfdata(i).y=rd{i}(:,:);   
    sz=size(mfdata(i).y,1);
    mfdata(i).y=(mfdata(i).y(1:end/2,:) + mfdata(i).y((end/2+1):end,:))/2;
    rd{i} = mfdata(i).y;    
    pts{i}=pts{i}(pts{i} <= sz/2);
    mfdata(i).y=mfdata(i).y(pts{i},:);
    mfdata(i).y=mfdata(i).y(:);
    mfdata(i).vy=(srd{i}(:,:)).^2;
    mfdata(i).vy=(mfdata(i).vy(1:end/2,:) + mfdata(i).vy((end/2+1):end,:))/4;
    mfdata(i).vy=mfdata(i).vy(pts{i},:);
    mfdata(i).vy=mfdata(i).vy(:);     
    mfdata(i).x=1:sz/2;
  else
    mfdata(i).y=rd{i}(pts{i},:);
    mfdata(i).y=mfdata(i).y(:); 
    mfdata(i).vy=(srd{i}(pts{i},:)).^2;
    mfdata(i).vy=mfdata(i).vy(:);
  end
  if 0
    modelfunc='gauss';
    modelguess=[1 .5 2e-3];
    modelmask=[1 1 1];
  else
    modelfunc='quad';
    modelguess=[1 -1e-3];
    modelmask=[1 1e-3];
  end
  
  [mfmodel(i).yfn mfmodel(i).fn] = skew_fit_function_polar([(12:14)+length(modelmask)*(i-1)],modelfunc);
end
    
% Add an extra fit group that forces singlet to line up.
if ~isempty(config.initial)
    p0 = config.initial;
else
  if minus    
    error('Please don''t do this');
    p0 =   [0, 0, 1, pi/2, 0, 1, .45*pi, pi/2, 1, -1, 1, repmat(modelguess,1,length(mfdata))];
  else    
    p0 =   [0, 0, 1, pi/2, pi, 1, pi/2 , pi/2, 1, -1, 1, repmat(modelguess,1,length(mfdata))];
  end
end
if ~isopt(config,'nofit')
    pmask= [0, 0, 0, 1   , 1 , 0, 1    ,  0  , 0 , 1, 0, repmat(modelmask,1,length(mfdata))];
    %
    [p, chisq] = mfitwrap(mfdata,mfmodel,p0,'plfit plotiter err',[0 0 0 0 0 0 0 0 0 0 0 repmat(modelmask,1,length(mfdata))] & pmask);
    fprintf('After decay fit, reduced chi^2 = %g\n',chisq);
    [p, chisq] = mfitwrap(mfdata,mfmodel,p,'plinit plfit plotiter err',[0 0 0 1 1 1 1 1 1 1 1 repmat(modelmask*0,1,length(mfdata))] & pmask);
    fprintf('After initial axis fit, reduced chi^2 = %g\n',chisq);
    [p, chisq] = mfitwrap(mfdata,mfmodel,p ,'plfit plotiter err',[0 0 0 1 1 1 1 1 1 1 1 repmat(modelmask,1,length(mfdata))] & pmask);
    fprintf('After overall fit, reduced chi^2 = %g\n',chisq);
    
    if 1 % fit losses
        pmask(6)=1;
        pmask(9)=1;
        [p, chisq] = mfitwrap(mfdata,mfmodel,p ,'robust plfit plotiter err', pmask);
        [p, chisq, p_cov] = mfitwrap(mfdata,mfmodel,p ,'plfit plotiter fine err', pmask);
        fprintf('After loss fit, reduced chi^2 = %g\n',chisq);
    end
    
    %[p, chisq] = mfitwrap(mfdata,mfmodel,p ,'plfit err plotiter fine', pmask);
    %fprintf('After final non-robust fit, reduced chi^2 = %g\n',chisq);
else
  p=p0;
end
  
pm=skew_fit_function_polar(p,'basis');
pm=reshape(pm(1:9),3,3);
pmi=inv(pm);
po=p(10:11);
% sensor correction
%a=-s/2 + t/2;
%b=s/2 + t/2;
%yp = y * a + b;
data_corrector=make_data_corrector(po,pm);
%
mfile=[rtfile(1:(end-4)) '_mat.mat'];
if ~isopt(config,'nofit')
  save(mfile,'pm','po','p_cov','data_corrector');
end


% Set up a powerpoint
ppt=guidata(pptplot);     
set(ppt.e_file,'String',rtfile);   
set(ppt.e_figures,'String','[60 61 255 256]');
set(ppt.e_title,'String',regexprep(regexprep(rtfile,'(sm_)|(\.mat)',''),'.*\',''));
if isopt(config,'nofit')
  set(ppt.e_body,'String','');     
else    
  bs=['M: ' sprintf('%+6.3f %+6.3f %+6.3f\n',pm) sprintf('\nP_0: ') num2str(po) sprintf('\np: '),num2str(p(1:11))];
  bs=[bs sprintf('\nP0: '),num2str(p0(1:11)),sprintf('\n Pmask'),num2str(pmask)];
  set(ppt.e_body,'String',bs);     
end

  

% dbloch style
%
c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','sphere','color',[c .05])];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|S>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|T>')];
cmds=[cmds struct('type','label','val',[0 0 1.15],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|\downarrow\uparrow>')];
if 1 % Add fit axes
  cmds=[cmds struct('type','vector','val',pmi(:,1),'size',1,'color',[0.0 0.0 0.0],'label','X')];
  cmds=[cmds struct('type','vector','val',pmi(:,2),'size',1,'color',[0.0 0.0 0.0],'label','Y')];
  cmds=[cmds struct('type','vector','val',pmi(:,3),'size',1,'color',[0.0 0.0 0.0],'label','Z')];  
end
opts=struct('opts',' resize raise','dt',.1);
base=cmds;
opengl software;


if 1 && ~strcmp(computer,'MACI64') % povray
c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
pcmds={};
pcmds=[pcmds struct('type','disc','val',[1 0 0],'color',[c .05],'color2',c)];
pcmds=[pcmds struct('type','disc','val',[0 1 0],'color',[c .05],'color2',c)];
pcmds=[pcmds struct('type','disc','val',[0 0 1],'color',[c .05],'color2',c)];
pcmds=[pcmds struct('type','sphere','color',[c .01])];
bra=myunicode('007c');
ket=myunicode('232a');
pcmds=[pcmds struct('type','label','val',[1.15 0 0],'label',[bra 'S' ket])];
pcmds=[pcmds struct('type','label','val',[-1.15 0 0],'label',[bra 'T' ket])];
pcmds=[pcmds struct('type','label','val',[0 0 1.15],'label',[bra myunicode('2191') myunicode('2193') ket])];
pcmds=[pcmds struct('type','label','val',[0 0 -1.15],'label',[bra myunicode('2193') myunicode('2191') ket])];
pcmds=[pcmds struct('type','vector','val',pmi(:,1),'size',1,'color',[1.0 0.0 0.0])];
pcmds=[pcmds struct('type','vector','val',-pmi(:,2),'size',1,'color',[1.0 0.0 0.0])];
pcmds=[pcmds struct('type','vector','val',pmi(:,3),'size',1,'color',[1.0 0.0 0.0])];  
popts=struct('dt',.1,'pbloch',1,'view',[.6, 6, .4]);
figure(255);
clf; 
pbloch(pcmds,popts);
else
figure(255);
clf; 
dbloch(cmds,opts);
end

figure(256);
clf; 
dbloch(cmds,opts);
hold on; 
samples=1000;
colors='rgbcmkk';
if 1
  markersize=5;
  markerline=1;
  linewidth=.5;
  errorlinewidth=2;
else
  markersize=20;
  markerline=2;
  linewidth=2;
  errorlinewidth=4;
end    

% datasets here.
datasets=1:length(rd)
%datasets=2
for i=datasets
    rdp{i}=rd{i};
    if 1
      rdp{i}=data_corrector(rd{i});
    else
      rdp{i}(:,2)=-rdp{i}(:,2);
    end
    %rdp{i}(:,1)=-rdp{i}(:,1);
    %rdp{i}(:,3)=-rdp{i}(:,3);
    
    %rdp{i}=rd{i}*pm;
    %plot3(rd{i}(:,1),rd{i}(:,2),rd{i}(:,3),[colors(i) '-'])
    f=plot3(rdp{i}(1,1),rdp{i}(1,2),rdp{i}(1,3),['kx']);
    set(f,'MarkerSize',markersize*5);
    set(f,'LineWidth',markerline*2);
    f=plot3(rdp{i}(1,1),rdp{i}(1,2),rdp{i}(1,3),[colors(i) 'x']);
    set(f,'MarkerSize',markersize*5);
    set(f,'LineWidth',markerline);

    
    switch 1
        case 1 % no error bars
          f=plot3(rdp{i}(:,1),rdp{i}(:,2),rdp{i}(:,3),[colors(i) '.']);
          set(f,'MarkerSize',markersize);
        case 3 % Plot really slow error bars
          f=plot3([(rdp{i}(:,1)-srd{i}(:,1)) (rdp{i}(:,1)+srd{i}(:,1))]' ...
            ,[rdp{i}(:,2) rdp{i}(:,2)]',[rdp{i}(:,3) rdp{i}(:,3)]',[colors(i) '-']); set(f,'LineWidth',errorlinewidth);
          f=plot3([(rdp{i}(:,1)) (rdp{i}(:,1))]' ...
            ,[rdp{i}(:,2)-srd{i}(:,2) rdp{i}(:,2)+srd{i}(:,2)]',[rdp{i}(:,3) rdp{i}(:,3)]',[colors(i) '-']); set(f,'LineWidth',errorlinewidth);
          f=plot3([(rdp{i}(:,1)) (rdp{i}(:,1))]' ...
            ,[rdp{i}(:,2) rdp{i}(:,2)]',[rdp{i}(:,3)-srd{i}(:,3) rdp{i}(:,3)+srd{i}(:,3)]',[colors(i) '-']); set(f,'LineWidth',errorlinewidth);
        case 3 % plot error ellipsoids
          ecmds={};
          for j=1:size(covrd,2)
             ecmds{end+1}=struct('type','ellipse','sigma',covrd{i,j},'center',rdp{i}(j,:),'color',[1 0 0 .3]);
          end
          eopts=opts;         
          dbloch(ecmds,struct('opts','hold'));
    end
    
    %intf=@pchip; % interpolation function
    intf=@spline; % interpolation function
    switch 2
        case 1 % nothing
        case 2 % splines
            intt=linspace(1,size(rd{i},1),samples);
            t=1:size(rd{i},1);
            f=plot3(intf(t,rdp{i}(:,1),intt),intf(t,rdp{i}(:,2),intt),intf(t,rdp{i}(:,3),intt),[colors(i) '-']);
            set(f,'LineWidth',linewidth);
        case 3
            f=plot3(rdp{i}(:,1),rdp{i}(:,2),rdp{i}(:,3),[colors(i) '-']);
            set(f,'LineWidth',linewidth);
    end
end
%set(gca,'xlim',[-1.1 1.1]);
%set(gca,'ylim',[-1.1 1.1]);
%set(gca,'zlim',[-1.1 1.1]);
xlabel(axismap{1});
ylabel(axismap{2});
zlabel(axismap{3});
%l=legend(cellfun(@num2str,num2cell(datasets),'UniformOutput',false))
%l=legend(rt.scan.data.pulsegroups(r1s(1)).name,rt.scan.data.pulsegroups(r2s(1)).name)
%set(l,'interpreter','none');
return;

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;

function data_corrector=make_data_corrector(po,pm)
data_corrector=@(data) (data*(-po(1)/2+po(2)/2)+po(1)/2+po(2)/2)*pm;
return;