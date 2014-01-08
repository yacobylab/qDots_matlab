%% Try to unscrambles axes based on experimental data.
% sm_TomoCal_rot_R_5281.mat is nice 3-J axis
% Load the trivial tomography
tt=0;
if tt
if ~exist('ttfile') | ttfile == 0
    ttfile='a';
end
end
if ~exist('rtfile') | rtfile == 0
    rtfile='a';
end
if tt
  ttfile=uigetfile('sm_TomoCal_trivial*.mat','Trivial Set',ttfile);
end
[rtfile dir]=uigetfile('sm_TomoCal_rot*.mat','Rotating Set',rtfile);
rtfile = [dir rtfile];
if tt
  tt=load(ttfile);
end
rt=load(rtfile);
scantime=getscantime(rt.scan,rt.data);

% Make some color plots for switch checks
figure(1);
clf;
d=permute(rt.data{1},[1 3 2]);
imagesc(reshape(d,size(d,1),size(d,2)*size(d,3)));
%% Scaling for non-raw calibration data
if ~isempty(strfind(rtfile,'_L_'))
    side='left';
else
    side='right';
end
%[t1t t1] = att1(side,scantime,'before');
t1=0.05;
%[t1t t1] = att1(side,scantime,'after');
%t1 = 4.5/55;  % t1 on left.
[rtd rtff mv fp]=anaHistScale(rt.scan,rt.data,t1,2:length(rt.scan.data.pulsegroups));
rtd=rtd{1}*2-1;
rtff=rtff{1};
if tt
  ttd=cellfun(rtff,tt.data,'UniformOutput',false);
  ttd=ttd{1}*2-1;
end

% Definitions to help us keep axes straight. 

axismap={'ST','Y','UD'};
axisorder=[1 3 2];
 
%% Figure out the angle to the singlet-triplet axis using the trivial tomography.
% find the null-rotation in the TT group
if tt
  nulls=find(~cellfun(@(x)isempty(x),strfind({tt.scan.data.pulsegroups.name},'null')));
  % pull out the appropriate data
  nvec=nanmean(ttd(:,nulls(axisorder)),1)';
  nvecs=nanstd(ttd(:,nulls(axisorder)),1,1)'/sqrt(size(ttd,1));
  ybasis=[nvec(2)/nvec(1) 1 0];
  udbasis=[nvec(3)/nvec(1) 0 1];
  stbasis=[1 0 0];
  basis0=[stbasis' ybasis' udbasis'];
end

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
  %pts{i}=1:32;
  pts=1:length(mfdata(i));
  mfdata(i).x=pts{1};

  xvals=plsinfo('xval',rt.scan.data.pulsegroups(rs{i}(1)).name);
  if xvals(end/2+1) == xvals(1)  % This is folded data...
    fprintf('Group %d appears to be folded\n',i);
    mfdata(i).y=rd{i}(:,:);   
    sz=size(mfdata(i).y,1);
    mfdata(i).y=(mfdata(i).y(1:end/2,:) + mfdata(i).y((end/2+1):end,:))/2;
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
if minus    
  p0 =   [0, 0, 1, pi/2, 0, 1, .45*pi, pi/2, 1, -1, 1, repmat(modelguess,1,length(mfdata))];
else
  p0 =   [0, 0, 1, pi/2, pi, 1, pi/2 , pi/2, 1, -1, 1, repmat(modelguess,1,length(mfdata))];
end
  pmask= [0, 0, 0, 1   , 1 , 0, 1    ,  1  , 1 , 1 ,0, repmat(modelmask,1,length(mfdata))];
%
[p, chisq] = mfitwrap(mfdata,mfmodel,p0,'plfit err plotiter',[0 0 0 0 0 0 0 0 0 0 0 repmat(modelmask,1,length(mfdata))] & pmask);
fprintf('After decay fit, reduced chi^2 = %g\n',chisq);
[p, chisq] = mfitwrap(mfdata,mfmodel,p,'plinit plfit err plotiter',[0 0 0 1 1 1 1 1 1 1 1 repmat(modelmask*0,1,length(mfdata))] & pmask);
fprintf('After initial axis fit, reduced chi^2 = %g\n',chisq);
[p, chisq] = mfitwrap(mfdata,mfmodel,p ,'plfit err plotiter',[0 0 0 1 1 1 1 1 1 1 1 repmat(modelmask,1,length(mfdata))] & pmask);
fprintf('After overall fit, reduced chi^2 = %g\n',chisq);

if 1 % fit losses
  pmask(6)=1;
  pmask(9)=1;
  [p, chisq] = mfitwrap(mfdata,mfmodel,p ,'robust plfit err plotiter', pmask);
  [p, chisq] = mfitwrap(mfdata,mfmodel,p ,'plfit err plotiter fine', pmask);
  fprintf('After loss fit, reduced chi^2 = %g\n',chisq);
end

%[p, chisq] = mfitwrap(mfdata,mfmodel,p ,'plfit err plotiter fine', pmask);
%fprintf('After final non-robust fit, reduced chi^2 = %g\n',chisq);
  
pm=skew_fit_function_polar(p,'basis');
pm=reshape(pm(1:9),3,3);
pmi=inv(pm);
po=p(10:11);
% sensor correction
%a=-s/2 + t/2;
%b=s/2 + t/2;
%yp = y * a + b;
data_corrector=@(data) (data*(-po(1)/2+po(2)/2)+po(1)/2+po(2)/2)*pm;
pm
po
%%
mfile=[rtfile(1:(end-4)) '_mat.mat'];
save(mfile,'pm','po','data_corrector');


%% Set up a powerpoint
ppt=guidata(pptplot);     
set(ppt.e_file,'String',rtfile);   
set(ppt.e_figures,'String','[60 61]');
set(ppt.e_title,'String','Tomography Calibration');
bs=['M: ' sprintf('%+6.3f %+6.3f %+6.3f\n',pm) sprintf('\nP_0: ') num2str(po) sprintf('\np: '),num2str(p(1:11))];
bs=[bs sprintf('\nP0: '),num2str(p0(1:11)),sprintf('\n Pmask'),num2str(pmask)];
set(ppt.e_body,'String',bs);     

  

%% oversampled plot, spline interpolation
%% dbloch style
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
scale = 5;
dbz=[0 0 1/32];
j1 =[1/9 0 0]+dbz;
atan2(j1(3),j1(1))
j2 =[norm(dbz) 0 0]+dbz;
dbzp=struct('type','vector','val',dbz*scale,'size',1,'color',[0.5 0.5 0.1],'label','\delta B_z');
jp0=struct('type','vector','val',(j1-dbz)*scale,'size',1,'color',[1 0 0],'label',' J');
jp=struct('type','vector','val',j1*scale,'size',1,'color',[1 0 0],'label',' Rotation Axis');
jpp=struct('type','vector','val',j2*scale,'size',1,'color',[0 0 1],'label',' J''');

opts=struct('opts',' resize raise','dt',.1);
base=cmds;
opengl software;


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

%% oversampled plot, spline interpolation
%% pbloch style

%% set up sphere
c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
c2=1-(1-c)/5;
c1=c;
c=1-(1-c)/2;
%c=[200/255 255/255 75/255];
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c2 .001],'color2',[c1 0.01],'trans','rgbf')];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c2 .001],'color2',[c1 0.01],'trans','rgbf')];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c2 .001],'color2',[c1 0.01],'trans','rgbf')];
cmds=[cmds struct('type','sphere','color',[c .005],'trans','rgbf')];
bra=myunicode('007c');
ket=myunicode('232a');
%cmds=[cmds struct('type','label','val',[1.15 0 0],'label',[bra 'S' ket])];
%cmds=[cmds struct('type','label','val',[-1.15 0 0],'label',[bra 'T' ket])];
%cmds=[cmds struct('type','label','val',[0 0 1.15],'label',[bra myunicode('2191') myunicode('2193') ket])];
%cmds=[cmds struct('type','label','val',[0 0 -1.15],'label',[bra myunicode('2193') myunicode('2191') ket])];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label',[bra 'Z' ket])];
%cmds=[cmds struct('type','label','val',[-1.15 0 0],'label',[bra 'T' ket])];
cmds=[cmds struct('type','label','val',[0 0 1.15],'label',[bra 'X' ket])];
%cmds=[cmds struct('type','label','val',[0 0 -1.15],'label',[bra myunicode('2193') myunicode('2191') ket])];

scale = 5;
opts=struct('dt',.1,'pbloch',1,'view',[.6, 6, .4]);

colors=[1 0 0 ; .1 .5 0 ; 0 0 1 ; 0 1 1 ; 1 0 1; .5 .5 0; 0 0 0];

if 1
  markersize=.0125;
  linewidth=.0045;
  errorlinewidth=2;
end
tcmds=cmds;
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
    tcmds=[tcmds struct('type','sphere','radius',markersize,'value',rdp{i},'color',colors(i,:))];
    tcmds=[tcmds struct('type','spline','radius',linewidth, 'value',rdp{i},'color',[colors(i,:) .5])];
end

pbloch(tcmds,opts)
%% As above, but render a rotation
mv={};
opts=struct('dt',.1,'pbloch',1,'view',[.6, 6, .4]);
r=norm(opts.view);
phi=atan2(opts.view(3),opts.view(2));
theta=atan2(opts.view(1), norm(opts.view(2:3)));
m=[];
sequence.framemult=18;
p=linspace(phi,phi+2*pi,10*sequence.framemult);
t=linspace(theta,theta,10*sequence.framemult);
for i=1:length(p)
    opts.view=r*[sin(t(i)), cos(t(i)) * cos(p(i)), cos(t(i)) * sin(p(i))];   
    %opts.right=[0, cos(p(i)-pi/2), sin(p(i)-pi/2)];   
    %opts.up=-cross(opts.view,opts.right);
if 0
    for j=1:length(tcmds)
           if isfield(tcmds{j},'val') && (norm(cmds{j}.val-[0 0 -1.15]) < 1e-6 || norm(cmds{j}.val-[-1.15 0 0]) < 1e-6)
               tcmds{j}.color(4)=1-dot(opts.view,cmds{j}.val)^2;
           end
    end               
end
    [po img] = pbloch([tcmds],opts);
    m = [ m img ] ;
end
mv = [ mv m];

