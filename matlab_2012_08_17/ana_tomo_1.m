function out = ana_tomo_1(fname, config)
%function out = ana_tomo_1(fname, config)

% analyze one-qubit tomography data.  
% can handle scans with tomography of more than one rotation. 
% assumes pulsegroups have the names of the tomo (ST, UD, Y)
% guesses whether groups are pack with tomo changing on inside loop or
% outside loop.
% data will come out in X-Y-Z order, along wiht 

if ~exist('fname','var') || isempty(fname)
   fname = uigetfile('sm*tomo*.mat'); 
end

if isempty(fname)|| (isnumeric(fname)) % hack because empty file gives 0
    warning('no file found \n');
    return
end

if ~exist('config','var')
    config=struct();
end
% set up default config stuff
config = def(config, 'opts','pauli');   % Random boolean options
config = def(config, 'figbase',700); 
config = def(config, 'Ymult',1); 
config = def(config, 'calFile',''); 

side = getside(fname);
s=load(fname); %load the data
plsgrps = {s.scan.data.pulsegroups.name};
scantime=getscantime(s.scan,s.data);
[t1t t1] = att1(side,scantime,'ask' );
[axisorder, dsets, names, dbzgrps,tomogrps]=ordertomo(plsgrps); %extract all the groups
%the above should come out in order ST,UD,Y
config = def(config, 'groups',tomogrps); 
[rtd rtff mv fp]=anaHistScale(s.scan,s.data,t1,config.groups); % scale the data
rtd= -1*(rtd{1}*2-1);
rtff=rtff{1};
d = squeeze(nanmean(rtd,1));

%now time calibrate tomography
if isempty(config.calFile)
    str = ['sm*_', upper(side(1)),'*_mat*'];
    fprintf('find tomoCal file \n');
    calfile = uigetfile(str);
    calname = calfile;
    if ~ischar(calfile) %protect against empty files: defaults to no calibration
        %calfile.pm = eye(3);
        %calfile.po = [-1 1];
        %calfile.data_corrector= @(data)(data*(-po(1)/2+po(2)/2)+po(1)/2+po(2)/2)*pm;
        calfile.data_corrector = @(x) x;
        warning('no tomocal file. assuming no correction');
    else
        calfile = load(calfile);
    end
else
    calfile = load(config.calFile{1});
end
% apply sensor correction
%sensorcorrection = @(x,po) po(1)+.5*(po(2)-po(1))*(x+1);



rdp = cell(1,length(dsets));
permuter = [2 3 1];% change normal ST,UD,Y into X,Y,Z
for j = 1:length(dsets)
   rdp{j} =  d(dsets{j}(permuter),:);
   %rdp{j}=rdp{j}(:,permuter);
    if ~isopt(config, 'nosensorcorrection')
      %rdp{j} = sensorcorrection(rdp{j},calfile.po)'*calfile.pm;
      %good-data = pm * messed-up-data
      %rdp{j} = (calfile.pm*sensorcorrection(rdp{j},calfile.po))';
      rdp{j} = calfile.data_corrector(rdp{j}');
    end
    
end
out.data = rdp;
out.plsgrps = plsgrps;
out.dsets = dsets;
out.axisorder =axisorder;
tomos = {'ST','UD','Y'};
out.tomos = tomos(permuter);
out.name = fname;
out.cal = calname;

if isopt(config,'pauli')
    for j = 1:length(rdp)
        xvals = 1:length(rdp{j});
        figure(config.figbase);
        config.figbase =config.figbase+1;
        clf; hold on;
        l=size(rdp{j},2);
        for i=1:l
            area(xvals, rdp{j}(:,i)+2*(l-i), 2*(l-i),'FaceColor',[1 0 0]);
        end
        set(gca,'YTick',2*(0:1:(l-1)));
        set(gca,'YTickLabel',out.tomos(end:-1:1));
        set(gca,'YLim',[-1 l*2-1]);
        set(gcf,'Name','Pauli Plot');
    end
end

%now complicated plotting
if ~isopt(config,'noplot')
    axismap = tomos(axisorder(permuter));
    figs = [];
    c=[255/255 139/255 29/255];
    c=[139/255 255/255 29/255];
    [cmds, opts, base]=makecmds(c);
    opengl software;
    
%     figure(config.figbase);
%     figs=[figs config.figbase]; config.figbase=config.figbase+1;
%     clf;
%     dbloch(cmds,opts);
%     
%     
    figure(config.figbase);
    figs=[figs config.figbase]; config.figbase=config.figbase+1;
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
    datasets=1:length(dsets);
    badinds = {};
    for i=datasets
        
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
        
        %plot big markers on points outside of Bloch Sphere
        norms = sum(rdp{i}'.^2); %norms;
        badinds{i} = find(norms>1.05); 
        f=plot3(rdp{i}(badinds{i},1),rdp{i}(badinds{i},2),rdp{i}(badinds{i},3),'k^');
        set(f,'MarkerSize',markersize*3);
        set(f,'LineWidth',markerline*2);
        
        %intf=@pchip; % interpolation function
        intf=@spline; % interpolation function
        switch 2
            case 1 % nothing
            case 2 % splines
                intt=linspace(1,size(rdp{i},1),samples);
                t=1:size(rdp{i},1);
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
    %legend(names); %doesnt work
    %set(l,'interpreter','none');
    out.badinds = badinds;
end

% Set up a powerpoint
if ~isopt(config,'noppt')
ppt=guidata(pptplot);     
set(ppt.e_file,'String',fname);   
set(ppt.e_figures,'String',['[' sprintf('%d ',figs) ']']);
set(ppt.e_title,'String',regexprep(regexprep(fname,'(sm_)|(\.mat)',''),'.*\',''));
descr = sprintf('single qubit tomography \n datasets \n');
descr = [descr, sprintf('%s \n',names{:})];
set(ppt.e_body,'String',descr);
clipboard('copy',descr);
end

end

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

function side = getside(name)
%sort out the side
isright = (strfind(name,'_R_'));
isleft = (strfind(name,'_L_'));
if isempty(isright) && isempty(isleft)
    error('cannot determine side');
elseif isempty(isright)
    side = 'left';
elseif isempty(isleft)
    side = 'right';
else %found _R and _L, default to the first one
    if isright > isleft
        side = 'left';
    else
        side = 'right';
    end
end
end

function [axisorder, dsets, names, dbzgrps, tomogrps]=ordertomo(plsgrps)
tomos = {'ST','UD','Y'};
dbzgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps,'^([dD][bB][zZ])')));
if length(dbzgrps)==length(plsgrps) %all the groups are dbz
   tomogrps = dbzgrps; 
else
  refgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps,'^(ref)'))); %deprecated, hopefully empty
  tomogrps = setdiff(1:length(plsgrps),union(dbzgrps,refgrps));
end
%some sanity checks
switch mod(length(plsgrps),3)
    case 0
        fprintf('found multiple of 3 groups \n');
    case 1
        if ~isempty(dbzgrps)
            fprintf('found 3n+dbz groups \n');
        else
            error('found 3n+1 groups but non with dbz reference');
        end
    otherwise
        error('Number of groups is 3 * integer +2. Confusing... \n')
end

%STgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps(tomogrps),'^([sS][tT])')));
STgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps(tomogrps),'ST')));
UDgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps(tomogrps),'UD')));
Ygrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps(tomogrps),'Y')));

%get the axis order. Assumed to be the same.
[dummy axisorder] = sort([STgrps(1), UDgrps(1), Ygrps(1)]);

if length(STgrps) ==1 %just one set of tomography
    dsets = {permute(tomogrps,axisorder)};
    ii = strfind(plsgrps{STgrps},'ST');
    names = {plsgrps{STgrps}([1:ii-1, ii+2:end])};
else %more than one tomography group
    
    %find whether groups are packed with tomo changing on inner loop or outter
    %loop
    switch diff(STgrps(1:2));
        case 1
            innerTomo = 0;
        case 3
            innerTomo = 1;
        otherwise
            error('cannot determine packing order of groups');
    end
    dsets = {}; names = {};
    if innerTomo
        for j = 1:3:length(tomogrps)
            dsets{end+1} =  tomogrps(j)+permute((0:2),axisorder); %ggrp(1)+(j:j+2)-1;
            ii = strfind(plsgrps{tomogrps(j)},tomos{axisorder(1)});
            names{end+1} = plsgrps{tomogrps(j)}([1:ii-1, ii+2:end]);
        end
    else %tomography is swept in outter loop
        for j = 1:length(tomogrps)
            dsets{end+1} = tomogrps(permute(j:3:end,axisorder));
            ii = strfind(plsgrps{tomogrps(j)},tomos{axisorder(1)});
            names{end+1} = plsgrps{tomogrps(j)}([1:ii-1, ii+2:end]);
        end
        
    end
end
end

function [cmds, opts, base]=makecmds(c)
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','sphere','color',[c .05])];
cmds=[cmds struct('type','label','val',[0 0 1.15],'label','|S>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|T>')];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|\downarrow\uparrow>')];
opts=struct('opts',' resize raise','dt',.1);
base=cmds;
end