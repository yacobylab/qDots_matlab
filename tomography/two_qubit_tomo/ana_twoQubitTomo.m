function val = ana_twoQubitTomo(file, config)
%function val = ana_twoQubitTomo(file, config)
%analyze two qubit tomo data
%file is a filename, left empty it prompts uigetfile
% multiselect is on. for now multiple files will just concatinate the data
%config is a struct with the usual format
%config.t1Time is a cell array that dictates whether T1 data is taken from
%before or after. elements can be 'before' or 'after'. defaults to 
%{'before', 'before');
% config.Ymult is a multiplier for the Y data- useful if not using tomocal
% config.opts can include 'fid' for bell state fidelity, 'angles' 
% to get the correct bell state angles. 
% opts can have 'noppt' to not fill out pptplot
% the function will prompt the user to load tomocal files (not the data
% files but the results of tomocal (pm, po, etc) or to ignore it

global tunedata;
%get file if needed
if ~exist('file','var') || isempty(file)
    file = uigetfile('sm*tomo*.mat', 'MultiSelect', 'on');
end

if ~iscell(file)
    file = {file};
end

if ~exist('config','var')
    config=struct();
end
% set up default config stuff
config = def(config,'opts','fid');   % Random boolean options
config = def(config,'t1Time',{'before', 'before'}); % defaults to try to get t1 from before {left, right}
config = def(config, 'Ymult', []); % default Y-multiplier = 1; can set to -1 if not calibrating tomo. empty=auto


for j = 1:length(file)
fprintf('%s \n', file{j});
end

rd = {};

%load files
for j = 1:length(file)
s=load(file{j});
rdr=anaRawUnpack(s.scan,s.data);
scantime=getscantime(s.scan,s.data);



% histogram the data
%first get t1
[t1lt t1l] = att1('left',scantime,config.t1Time{1});
[t1rt t1r] = att1('right',scantime,config.t1Time{2});

%histogram it
rdd=anaRawScale(rdr,[t1l t1r],2:length(s.scan.data.pulsegroups));
for i=1:2
    rdd{i}=-(rdd{i}*2-1);
end

for k = 1:2; 
    if j==1
        rd{k} = rdd{k};
    else
    rd{k} = [rd{k}; rdd{k}]; 
    end
end

end
corrected =0;
figind = 0;
bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};

switch (input('Correct Tomography Errors (yes/[no])? ','s'))
    case 'yes'
        if isempty(config.Ymult)
            config.Ymult=1;
        end
        
        fprintf('find Left tomoCal file \n');
        ssL = load(uigetfile('sm*L*_mat*'));
        fprintf('find Right tomoCal file \n');
        ssR = load(uigetfile('sm*R*_mat*'));
        % apply sensor correction
        sensorcorrection = @(x,po) po(1)+.5*(po(2)-po(1))*(x+1);
        
        if ~isopt(config, 'nosensorcorrection')
            rd{1} = sensorcorrection(rd{1},ssL.po);
            rd{2} = sensorcorrection(rd{2},ssR.po);
        end
        
        
        % ST U/D Y in file.  Right sweeps faster.
        % ST/ST(2) ST/UD(3) ST/Y(4)
        % UD/ST(5) UD/UD(6) UD/Y(7)
        %  Y/ST(8)  Y/UD(9)  Y/Y(10)
        
        
        IX=  squeeze(nanmean(nanmean(rd{2}([2 5 8],:,:),3),1))';
        IY=  config.Ymult*squeeze(nanmean(nanmean(rd{2}([4 7 10],:,:),3),1))';
        IZ=  squeeze(nanmean(nanmean(rd{2}([3 6 9],:,:),3),1))';
        XI=  squeeze(nanmean(nanmean(rd{1}([2 3 4],:,:),3),1))';
        YI=  config.Ymult*squeeze(nanmean(nanmean(rd{1}([8 9 10],:,:),3),1))';
        ZI=  squeeze(nanmean(nanmean(rd{1}([5 6 7],:,:),3),1))';
        XX =              nanmean(squeeze(rd{1}(2,:,:).*rd{2}(2,:,:)),2);
        XY = config.Ymult*nanmean(squeeze(rd{1}(4,:,:).*rd{2}(4,:,:)),2);
        XZ =              nanmean(squeeze(rd{1}(3,:,:).*rd{2}(3,:,:)),2);
        YX = config.Ymult*nanmean(squeeze(rd{1}(8,:,:).*rd{2}(8,:,:)),2);
        YY =              nanmean(squeeze(rd{1}(10,:,:).*rd{2}(10,:,:)),2);
        YZ = config.Ymult*nanmean(squeeze(rd{1}(9,:,:).*rd{2}(9,:,:)),2);
        ZX =              nanmean(squeeze(rd{1}(5,:,:).*rd{2}(5,:,:)),2);
        ZY = config.Ymult*nanmean(squeeze(rd{1}(7,:,:).*rd{2}(7,:,:)),2);
        ZZ =              nanmean(squeeze(rd{1}(6,:,:).*rd{2}(6,:,:)),2);
        
        %left and right rotation matricies
        ML = ssL.pm;
        MR = ssR.pm;
        
        %       1  2  3  4  5  6  7   8  9  10  11  12  13  14  15
        DDD = [XI YI ZI IX IY IZ XX, XY,XZ, YX, YY, YZ, ZX, ZY, ZZ];
        newD=0*DDD;
        
        L = [ML' zeros(3,12)];
        R = [zeros(3,3) MR' zeros(3,9)];
        MM = kron(ML',MR');
        M = [L; R; [zeros(9,6) MM]];
        
        for i = 1:length(XX)
            newD(i,:)  = M*DDD(i,:)';
        end
        
%        bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};
        dataexp=newD(:,[1:6 8:10 12:14 7 11 15])'; % reshape back into order of bs
        corrected = 1;
        
    case 'no'
        if isempty(config.Ymult)
            config.Ymult=-1;
        end
        
        % ST U/D Y in file.  Right sweeps faster.
        % ST/ST(2) ST/UD(3) ST/Y(4)
        % UD/ST(5) UD/UD(6) UD/Y(7)
        %  Y/ST(8)  Y/UD(9)  Y/Y(10)
        
        
        IX=  squeeze(nanmean(nanmean(rd{2}([2 5 8],:,:),3),1))';
        IY=  config.Ymult*squeeze(nanmean(nanmean(rd{2}([4 7 10],:,:),3),1))';
        IZ=  squeeze(nanmean(nanmean(rd{2}([3 6 9],:,:),3),1))';
        XI=  squeeze(nanmean(nanmean(rd{1}([2 3 4],:,:),3),1))';
        YI=  config.Ymult*squeeze(nanmean(nanmean(rd{1}([8 9 10],:,:),3),1))';
        ZI=  squeeze(nanmean(nanmean(rd{1}([5 6 7],:,:),3),1))';
        XX =              nanmean(squeeze(rd{1}(2,:,:).*rd{2}(2,:,:)),2);
        XY = config.Ymult*nanmean(squeeze(rd{1}(4,:,:).*rd{2}(4,:,:)),2);
        XZ =              nanmean(squeeze(rd{1}(3,:,:).*rd{2}(3,:,:)),2);
        YX = config.Ymult*nanmean(squeeze(rd{1}(8,:,:).*rd{2}(8,:,:)),2);
        YY =              nanmean(squeeze(rd{1}(10,:,:).*rd{2}(10,:,:)),2);
        YZ = config.Ymult*nanmean(squeeze(rd{1}(9,:,:).*rd{2}(9,:,:)),2);
        ZX =              nanmean(squeeze(rd{1}(5,:,:).*rd{2}(5,:,:)),2);
        ZY = config.Ymult*nanmean(squeeze(rd{1}(7,:,:).*rd{2}(7,:,:)),2);
        ZZ =              nanmean(squeeze(rd{1}(6,:,:).*rd{2}(6,:,:)),2);
        
        estring=['[' sprintf('%s(j) ',bs{:}) '];'];
        dataexp = zeros(15,size(XX,1));
        for j = 1:size(XX,1)
            dataexp(:,j) = eval(estring);
            
        end
    otherwise
        error('You Suck at Matlab');        
end

badxval = false;
paramsL = [];
paramsR = [];
  for i = 2:length(s.scan.data.pulsegroups)
        gd=plsinfo('gd', s.scan.data.pulsegroups(i).name, [],scantime);
        pL = plsinfo('params', gd.pulses.groups{1}, [],scantime);
        pR = plsinfo('params', gd.pulses.groups{2}, [],scantime);
        if isempty(pL)||isempty(pR)
            fprintf('Warining; no parameters on group %s\n',s.scan.data.pulsegroups(i).name);
            fprintf('Probable logging error.  Throwing out data.\n');
            pL = [0 0 0 0]; pR = pL;
            badxval = true;
        end
        paramsL=[paramsL; pL];
        paramsR=[paramsR; pR];
  end
  epsL = paramsL(:,1);
  epsR = paramsR(:,1);
  if any(diff(epsL)~=0) %is epsilon changing
      epsStringL = 'Could not determine epsL\n';
  else
      epsStringL = sprintf('epsL = %.2f \n', epsL(1));
  end

  if any(diff(epsR)~=0)
      epsStringR = 'Could not determine epsR\n';
  else
      epsStringR = sprintf('epsR = %.2f \n', epsR(1));
  end
  
  [r c]=find(diff(paramsL,[],2) ~= 0);
  if ~any(any(diff(paramsL,2))~=0)
      xvals = plsinfo('xval', gd.pulses.groups{1}, [],scantime);
      if size(xvals,1)==2 % 2D xval is usually sign of 2D qubit rotation scan
          xv2 = xvals;
          xvals = 1:size(xvals,2);
          xlab = 'Single Qubit Rotation';
          %config.opts = [config.opts, ' 2Dfid']; % make a 3d bar plot
          
      end
      if length(xvals)> size(s.data{1},3)
          %terrible hack to get the correct xvals. 
          xvals = xvals(end-1,:); 
          xlab = 'T_{entangle} (\mus)';
      end
      
  else
      
  switch mode(c)
      case 1
          xvals = paramsL(:,1);
          xlab = 'Epsilon_{entangle}';
      case 2
          xvals = paramsL(:,1);
          xlab = 'T_{entangle}';
      case 3
          xvals = paramsL(:,1);
          xlab = 'Epsilon_{rotation}';
      case 4
          xvals = paramsL(:,1);
          xlab = 'T_{entangle}';
      otherwise
        badxval = true;
        xlab = 'No friggin clue';
          
  end
  end
  
  
%find the optimal Bell States for all of the states in dataexp

if badxval
    xvals = 1:length(norms);
end
bestangles= [];
beststates= [];
bestfids= [];
for j=1:size(dataexp,2)
    [ba,bs,bf,statefunc] = BellStateFinder(dataexp(:,j));
    bestangles(j,:) = ba;
    beststates(j,:) = bs;
    bestfids(j) = bf;
    norms(j)=sqrt((sum(dataexp(:,j).^2)+1)/4);
end


if isopt(config, 'fid')
    figure(200+figind); clf;
    plot(xvals,bestfids,'rx-',xvals,norms,'bx-');
    ll = line([xvals(1),xvals(end)], [.5, .5]);
    set(ll, 'Color', [0 1 0]); %make it green
    set(ll, 'LineStyle', ':');
    title(sprintf('State Norm and Bell State Fidelity \n Maximum Fidelity = %.3f',max(bestfids)));
    xlabel(xlab);
    ylabel('Fidelity');
    legend('Bell State Fidelity', 'State Norm',...
        sprintf('Product states strictly forbidden'));
    yy= get(gca,'Ylim');
    if(yy(1) > .3)
        set(gca, 'Ylim', [.2 yy(2)]);
    end
    figind = figind +1;
end


if isopt(config, '2Dfid')
    maxfids2 = reshape(bestfids, round(sqrt(length(xvals))),round(sqrt(length(xvals)))); 
    figure(200+figind); clf;
    h =bar3(maxfids2);
shading interp;
set(gca,'zlim',[0 1]);
grid on;
for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k');    
    set(h,'FaceAlpha',.5);
    colormap(copper);
end

xlabel('Right Qubit Rotation');
ylabel('Left Qubit Rotation');
zlabel('Bell State Fidelity');
title(sprintf('Fidelity is between %.2f and %.2f', max(max(maxfids2)), min(min(maxfids2))));
    figind = figind +1;
    
    figure(200+figind);
    imagesc(maxfids2); colorbar;
    xlabel('Right Qubit Rotation'); ylabel('Left Qubit Rotation');
    figind = figind+1;
end

if isopt(config, 'angles')
    figure(200+figind); clf;
    subplot(2,1,1); hold on;
    [x y z] = sph2cart(bestangles(:,3), bestangles(:,1),1);
    plot3(x,y,z, 'rx');
    title('Left Qubit Best Angle');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    rotate3d on;
    subplot(2,1,2); hold on;
    [x2 y2 z2] = sph2cart(bestangles(:,4), bestangles(:,2),1);
    plot3(x2,y2,z2, 'rx');
    title('Right Qubit Best Angle');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    rotate3d on;
    
end

if ~isopt(config, 'noppt')
    flstr = '';
    for k = 1:length(file)
    flstr = [flstr, sprintf('%s \n',file{k})];
    end
     figs = 199+(1:figind);
     ppt=guidata(pptplot);     
     set(ppt.e_file,'String',file{1});     
     set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
     set(ppt.e_title,'String',regexprep(regexprep(file{1},'(sm_)|(\.mat)',''),'.*\',''));
     %set(ppt.e_title,'String','ana_twoQubitTomo');
     if corrected
     set(ppt.e_body,'String',[sprintf('two qubit state tomography (corrected)\n') epsStringL, epsStringR, flstr]);  
     else
         set(ppt.e_body,'String',[sprintf('two qubit state tomography (not corrected)\n') epsStringL epsStringR flstr]);  
     end
end

val = file;
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