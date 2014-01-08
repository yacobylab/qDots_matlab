function val = ana_twoQubitTomo_2(file, config)
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
% opts; noppt - no ppt windows setup
%       fid   - plot bell state fidelity
%      norm   - plot overall state norm
%       sqn   - plot single qubit norms (|<XI>+<YI>+<ZI>|)
%  concurrence - plot concurrence
%  eigens     - plot eigenvalues of rho.
%   fullplot  - provide full pauli plot.
%      sqt1   - plot non-oscillating amplitude in single qubits (<XI>,<IX>)
%      angles - plot fit bell state angles
%    fitfid   - attempt to fit data
%     stonly  - only rotate bell states around ST axis.
%       cal   - use tomocal
%       ana   - reload old result, if possible.
%       rng   - xval rng to fit.  default to [-inf inf]
%     subxbar - remove ST component of pauli plots.  reana needed.
% decoherence options;
%      sqwd - single qubit white dephasing
%      dqwd - single qubit white dephasing, different on each qubit
%               config.mfit : norm, sqn (single qubit norm), fid
% the function will prompt the user to load tomocal files (not the data
% files but the results of tomocal (pm, po, etc) or to ignore it
% FIXME; this help is woefully incomplete.

%get file if needed
if ~exist('file','var') || isempty(file)
    file = uigetfile('sm*tomo*.mat', 'MultiSelect', 'on');
end

if ~iscell(file)
    file = {file};
end
file=sort(file);

if ~exist('config','var')
    config=struct();
end
% set up default config stuff
config = def(config, 'opts','fid norm fidfit sqwd cal ana sqn fullplot angles concurrence');   % Random boolean options
config = def(config, 't1Time',{'before', 'before'}); % defaults to try to get t1 from before {left, right}
config = def(config, 'Ymult', []); % default Y-multiplier = 1; can set to -1 if not calibrating tomo. empty=auto
config = def(config, 'mfit','fid norm sqn'); % Things to fit in mfit mode.
config = def(config, 'rng', [-inf inf]);  % Range
config = def(config, 'calFile', '');  
bodystr='';
if isopt(config,'stonly')
    stonly='stonly';
else
    stonly='';
end

f=regexprep(file{end},'(\.mat)','');
if length(file) == 1 && ~isempty(strfind(file{1},'result'))
    savefile=file{1};
else
  if isopt(config,'cal')    
    savefile=sprintf('%s_results_cal.mat',f);
  else
    savefile=sprintf('%s_results.mat',f);
  end
end
fprintf('Save file is %s\n',savefile);
bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};

if isopt(config,'ana') && exist(savefile,'file')
    co=config;
    bso=bs;
    load(savefile);
    config=co;
    bs=bso;
else
    [dataexp covariances scantime s]=ana_TomoUnpack(file,config);
    
    figind = 0;    
    
    % Guess th correct xvals for the data.
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
    if badxval
        xvals = 1:length(norms);
    end
    bestangles= [];
    beststates= [];
    bestfids= [];
    %find the optimal Bell States for all of the states in dataexp

    [bestangles, beststates, bestpauli, bestfids_variances, bestfids, norms, norms_variances, statefunc, bpv] = fitBellStates(dataexp,covariances,stonly);
    s.data=[];    
    save(savefile);
end
if isopt(config,'reana')
    bestangles= [];
    beststates= [];
    bestfids= [];
    %find the optimal Bell States for all of the states in dataexp
    if isopt(config,'subxbar')
       xbl=mean(dataexp(1,5:end));
       z=zeros(1, size(dataexp,2));
       dataexp=dataexp-[repmat(xbl,1,size(dataexp,2)); z ; z ;z ;z ;z ;dataexp(5,:)*xbl ;dataexp(6,:)*xbl ;z;z;z;z ;xbl*dataexp(4,:);z;z];
       xbr=mean(dataexp(4,5:end));
       dataexp=dataexp-[z ; z ; z ; repmat(xbr,1,size(dataexp,2)) ; z ; z ; z ; z ; dataexp(2,:)*xbr ; z ; dataexp(3,:)*xbl ; z; dataexp(1,:)*xbr;z;z];
       bodystr=[bodystr sprintf('Subtracted <XI>=%g, <IX>=%g\n',xbl,xbr)];
    end
    [bestangles, beststates, bestpauli, bestfids_variances, bestfids, norms, norms_variances, statefunc,bpv] = fitBellStates(dataexp,covariances,stonly);
end

pts=find(xvals < config.rng(2) & xvals > config.rng(1));
val.epsL = epsL;
val.epsR = epsR;
% Handle two different versions of the norm.
if isopt(config,'oldnorm')
    for j=1:size(dataexp,2)
        norms(j)=sqrt((sum(dataexp(:,j).^2)+1)/4);        
        vweight=(-bpv'*bpv/norms(j)^3)/4;
        vweight=vweight-diag(diag(vweight));
        norms_variances(j)=sum(((dataexp(:,j)/norms(j)).^2) .*diag(covariances{j}),1)/16 - sum(sum(vweight .* covariances{j}));
    end
else
    for j=1:size(dataexp,2)
        norms(j)=sqrt((sum(dataexp(:,j).^2))/3);        
        vweight=(-bpv'*bpv/norms(j)^3)/3;
        vweight=vweight-diag(diag(vweight));
        norms_variances(j)=sum(((dataexp(:,j)/norms(j)).^2) .*diag(covariances{j}),1)/9 - sum(sum(vweight .* covariances{j}));
    end
end

% Make the "main" plot.
figure(200+figind); clf; hold on; lp=[];
if isopt(config,'fidfit')    
    style='x:';
else
    style='x-';
end
if isopt(config, 'fid')
    lp(end+1)=plot(xvals,bestfids,['m' style],'DisplayName','Fidelity');    
    if isopt(config,'err');
       p=errorbar(xvals,bestfids, sqrt(bestfids_variances),'m');              
    end
    ll = line([xvals(1),xvals(end)], [.5, .5]);
    lp(end+1)=ll;
    set(ll, 'DisplayName','Product States Forbidden');
    set(ll, 'Color', [0 1 0]); %make it green
    set(ll, 'LineStyle', ':');            
end

if isopt(config,'sqn')     % single qubit norms
    lp(end+1)=plot(xvals,sqrt(sum(dataexp(1:3,:).^2,1)),['r' style],'DisplayName','Left Norm');
    lp(end+1)=plot(xvals,sqrt(sum(dataexp(4:6,:).^2,1)),['b' style],'DisplayName','Right Norm');
end

if isopt(config,'sqt1')   % single qubit ST-axis projections. Should be constant
   lp(end+1)=plot(xvals,dataexp(1,:),'r:','DisplayName','Left <ST>');
   lp(end+1)=plot(xvals,dataexp(4,:),'b:','DisplayName','Right <ST>');   
end
if isopt(config,'norm')
    lp(end+1)=plot(xvals,norms',['k' style],'DisplayName','Norm');    
    if isopt(config,'err')
       errorbar(xvals,norms, sqrt(norms_variances),'b');       
    end
end

% Set up captions, etc.
if isopt(config,'err')
    [mf mi] = max(bestfids);
    title(sprintf('State Norm and Bell State Fidelity \n Maximum Fidelity = %.3f \\pm %.3f',mf,sqrt(bestfids_variances(mi))));
else
    title(sprintf('State Norm and Bell State Fidelity \n Maximum Fidelity = %.3f',max(bestfids)));
end
xlabel(xlab);
ylabel('Fidelity');
legend(lp);
yy= get(gca,'Ylim');
if(yy(1) > .3)
    set(gca, 'Ylim', [.2 yy(2)]);
end

if isopt(config, 'fidfit')
    % Returns coefficent of: left dot pv, right dot pv, off-axis left,
    % off-axis right, left/right terms.
    ncoeff=3;
    %!! comments below on axes are incorrect!!%% BEWARE!
    % <X> = S-T axis
    pfunc=@(p,x) [cos((x(:)+p(2))*p(3)) * p(1), ...  % <XI> - like coefficent
        cos((x(:)+p(2))*p(3)) * p(1), ...  % <IX> - like coefficent
        repmat(p(1)^2,length(x),1),   ...  % <XX> - like coefficent
        sin((x(:)+p(2))*p(3)) * p(1), ... % <XY> - like coefficent
        sin((x(:)+p(2))*p(3)) * p(1)]; % <YX> - like coefficent
    [m mi]=max(bestfids);
    cguess=[.9 0 pi/(2*xvals(mi))];
    fmtstring='';
    fitdescr='';
    fitdescr='Fidelity Fit';
    if isopt(config,'dqwd') % different single-qubit-white-dephasing, always seems to be singular...
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+2)), exp(-x(:)*(1/p(n+1)+1/p(n+2))),exp(-x(:)/p(n+1)),exp(-x(:)/p(n+2)) ];
        cguess(end+1:end+2)=[.5 2];
        fitdescr=[fitdescr ' 1Q white'];
        fmtstring=[fmtstring '\nSingle Qubit T_2^l: %g T_2^r: %g'];
    elseif isopt(config,'sqwd') % identical single-qubit-white-dephasing
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+1)), exp(-2*x(:)/p(n+1)),exp(-x(:)/p(n+1)),exp(-x(:)/p(n+1)) ];
        cguess(end+1)=[1];
        fmtstring=[fmtstring '\nSingle Qubit T_2: %g'];
        fitdescr=[fitdescr ' ident 1Q white'];
    elseif isopt(config,'dqgd') % different single-qubit gaussian dephasing
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-(x(:)/p(n+1)).^2), exp(-(x(:)/p(n+2)).^2), exp(-(x(:)*(1/p(n+1))).^2-(x(:)/p(n+2)).^2),exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+2)).^2) ];
        cguess(end+1:end+2)=[.5 2];
        fitdescr=[fitdescr ' 1Q white'];
        fmtstring=[fmtstring '\nSingle Qubit T_2^l: %g T_2^r: %g'];
    end
    
    if isopt(config,'sqgd') % identical single-qubit 1/f dephasing
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-(x(:)/p(n+1)).^2), exp(-(x(:)/p(n+1)).^2), exp(-2*(x(:)/p(n+1)).^2), ...
            exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2) ];
        cguess(end+1)=[1];
        fmtstring=[fmtstring '\nSingle Qubit 1/f T_2^*: %g'];
        fitdescr=[fitdescr ' ident. 1Q 1/f'];
    end
    if isopt(config,'tqgd') % two-qubit gate 1/f-noise dephasing
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2), (exp(-(x(:)/p(n+1)).^2)+ones(length(x),1))/2, ...
            exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2) ];
        cguess(end+1)=[1];
        fmtstring=[fmtstring '\nTwo Qubit 1/f T_2^*: %g'];
        fitdescr=[fitdescr ' 2Q 1/f'];
    end
    if isopt(config,'tqwd') % two-qubit gate white-noise dephasing
        n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-(x(:)/p(n+1)).^1),exp(-(x(:)/p(n+1)).^1), (exp(-(x(:)/p(n+1)).^1)+ones(length(x),1))/2, ...
            exp(-(x(:)/p(n+1)).^1),exp(-(x(:)/p(n+1)).^1) ];
        cguess(end+1)=[1];
        fmtstring=[fmtstring '\nTwo Qubit white T_2: %g'];
        fitdescr=[fitdescr ' 2Q white'];
    end
    
    if isopt(config,'tqmd') % Self-consistent model of independent epsilon noise on both qubits; basically YY doesn't dephase.
     n=length(cguess);
        pfunc=@(p,x) pfunc(p,x) .* [exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2), ones(length(x),1), ...
            exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2) ];
        cguess(end+1)=[1];
        fmtstring=[fmtstring '\nTwo Qubit modelled T_2: %g'];
        fitdescr=[fitdescr ' 2Q modeled'];       
    end
    
    fidfunc=@(p,x) transpose(1+sum(abs(subsref(pfunc(p,x),struct('type','()','subs',{{':',3:5}}))),2))/4;
    if isopt(config,'oldnorm')
      ampfunc=@(p,x) transpose(sqrt((1+sum(pfunc(p,x).^2,2))/4)); 
    else
      ampfunc=@(p,x) transpose(sqrt((sum(pfunc(p,x).^2,2))/3)); 
    end
    
    lnfunc=@(p,x) transpose(sqrt(sum(subsref(pfunc(p,x),struct('type','()','subs',{{':',1}})).^2,2)));
    rnfunc=@(p,x) transpose(sqrt(sum(subsref(pfunc(p,x),struct('type','()','subs',{{':',2}})).^2,2)));
    
    tmp.opts=config.mfit;
    mfdata=[];
    mfmodel=[];
    if isopt(tmp,'fid')
        mfdata(end+1).x=xvals(pts);
        mfdata(end).y=bestfids(pts);
        mfdata(end).name = 'fid';
        mfmodel(end+1).fn=fidfunc;
    end
    if isopt(tmp,'norm')
        mfdata(end+1).x=xvals(pts);
        mfdata(end).y=norms(pts);
        mfdata(end).name = 'norm';
        mfmodel(end+1).fn=ampfunc;
    end
    if isopt(tmp,'sqn')
        mfdata(end+1).x=xvals(pts);
        mfdata(end).y=sqrt(sum(dataexp(1:3,pts).^2,1));
        mfdata(end).name = 'Lnorm';
        mfmodel(end+1).fn=lnfunc;
        mfdata(end+1).x=xvals(pts);
        mfdata(end).y=sqrt(sum(dataexp(4:6,pts).^2,1));
        mfdata(end).name = 'Rnorm';
        mfmodel(end+1).fn=rnfunc;
    end
    [p, chisq] = mfitwrap(mfdata,mfmodel,cguess,'plfit plinit');
    val.p = p;
    val.chisq = chisq;
    val.mfdata = mfdata;
    val.dataexp = dataexp;
    val.angles = bestangles;
    val.pfunc = pfunc;
    val.fidfunc = fidfunc;

    x=linspace(xvals(1),xvals(end),512);
    figure(200+figind);
    hold on;
    if isopt(config,'fid')
       plot(x,fidfunc(p,x),'m-');
    end
    if isopt(config,'norm')
        plot(x,ampfunc(p,x),'k-');
    end
    if isopt(config,'sqn')
        plot(x,lnfunc(p,x),'r-');
        plot(x,rnfunc(p,x),'b-');
    end
    bodystr=[bodystr fitdescr sprintf('\n') sprintf('Amplitude: %g, W_ent: %g, T_ent: %g, Phase: %g',p([1 3]),pi/(2*p(3)), p(2))];
    bodystr=[bodystr sprintf(fmtstring,p(4:end))];
    
    if isopt(config,'gauss')
        bodystr=[bodystr sprintf('T2*_ent = %g\n',p(6))];
    end
    figind = figind +1;
    figure(200+figind);
    clf; hold on;
    plot(xvals*p(3)/pi,bestfids,'rx');
    plot(x*p(3)/pi,fidfunc(p,x),'k-');
end
figind = figind +1;

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
    set(gcf,'Name','2D Fidelity Plot');
end



if isopt(config,'concurrence')  % Check physicality of density matrix
  figure(200+figind); clf;
  figind=figind+1;
  c=[];
  mi=[];
  for i=1:size(dataexp,2)
       den=pauli2density(dataexp(:,i));
       e=eig(den);
       c(i)=real(concurrence(den));
       mi(i)=min(e);
  end
   lp=[];
   lp(end+1)=plot(xvals,c,'kx','DisplayName','Concurrence');
   ind=find(mi<0);
   if ~isempty(ind)
     hold on;
     lp(end+1)=plot(xvals(ind),c(ind),'r+','DisplayName','Unphysical \rho');
     set(lp, 'MarkerSize',8);
     set(lp, 'LineWidth',2);
   end
   legend(lp);
   xlabel(xlab);
   ylabel('Concurrence');
    set(gcf,'Name','Concurrences');   
end


if isopt(config,'eigens')  % Check physicality of density matrix
  figure(200+figind); 
  figind=figind+1;
  c=[];
  mi=[];
  ae=[];
  for i=1:size(dataexp,2)
       den=Pauli2density(dataexp(:,i));
       e=eig(den);
       mi(i)=min(e);
       ae(i,:)=sort(e);
  end
   clf;
   plot(xvals,mi,'kx','DisplayName','Minimum Eigenvalue');
   hold on;
   for i=2:size(ae,2)
     plot(xvals,ae(:,i),'rx','DisplayName','Other Eigenvalue');
   end
   title('Eigenvalues');
   xlabel(xlab);
   set(gcf,'Name','Eigenvalues');      
end

if isopt(config,'ptt')  % Check physicality of density matrix
  figure(200+figind); 
  figind=figind+1;
  c=[];
  mi=[];
  ae=[];
  for i=1:size(dataexp,2)
       den=Pauli2density(dataexp(:,i));
       e=eig(ptranspose(den));
       mi(i)=min(e);
       ae(i,:)=sort(e);
  end
   clf;
   plot(xvals,mi,'kx','DisplayName','Minimum Eigenvalue');
   hold on;
   for i=2:size(ae,2)
     plot(xvals,ae(:,i),'rx','DisplayName','Other Eigenvalue');
   end
   title('Eigenvalues');
   xlabel(xlab);
   set(gcf,'Name','Partial Transpose Test');      
end

if isopt(config, 'angles')
    figure(200+figind); clf;
    figind=figind+1;
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
    set(gcf,'Name','Angles');          
end
if isopt(config,'fullplot')
    figure(200+figind); 
    figind=figind+1;
    clf; hold on;
    colors=[1 0 0; 0 0 1 ; .7 0 .7];
    cind=[1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 ];
    l=size(dataexp,1);
    for i=1:l
        area(xvals, dataexp(i,:)+2*(l-i), 2*(l-i),'FaceColor',colors(cind(i),:));
    end
    set(gca,'YTick',2*(0:1:(l-1)));    
    set(gca,'YTickLabel',bs(end:-1:1));
    set(gca,'YLim',[-1 l*2-1]);
    set(gcf,'Name','Pauli Plot');          
end
if ~isopt(config, 'noppt')
    figs = 199+(1:figind);
    ppt=guidata(pptplot);
    set(ppt.e_file,'String',file{1});
    set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
    set(ppt.e_title,'String',regexprep(regexprep(file{1},'(sm_)|(\.mat)',''),'.*\',''));
    %set(ppt.e_title,'String','ana_twoQubitTomo');
    if  isopt(config,'cal')
        set(ppt.e_body,'String',[sprintf('two qubit state tomography (corrected)\n') epsStringL, epsStringR,bodystr]);
    else
        set(ppt.e_body,'String',[sprintf('two qubit state tomography (not corrected)\n') epsStringL epsStringR,bodystr]);
    end
end

val.file = file;
end


function [bestangles, beststates, bestpauli, bestfids_variances, bestfids, norms, norms_variances, statefunc, bpv] = fitBellStates(dataexp,covariances,stonly)
  for j=1:size(dataexp,2)
        [ba,bss,bf,bpv,statefunc] = BellStateFinder(dataexp(:,j),stonly);
        bestangles(j,:) = ba;
        beststates(j,:) = bss;
        bestpauli(j,:)=bpv;
        bestfids_variances(j)=sum((bestpauli(j,:)'/4).^2 .* diag(covariances{j}),1);
        bestfids(j) = bf;
        norms(j)=sqrt((sum(dataexp(:,j).^2)+1)/4);        
        vweight=(-bpv'*bpv/norms(j)^3)/4;
        vweight=vweight-diag(diag(vweight));
        norms_variances(j)=sum(((dataexp(:,j)/norms(j)).^2) .*diag(covariances{j}),1)/16 - sum(sum(vweight .* covariances{j}));        
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


function v=concurrence(m)
  sy=[0 -1i; 1i 0];
  m1=sqrtm(m);
  m2=kron(sy,sy);
  m3=sqrtm(m1*m2*(m.')*m2*m1);
  e=sort(eig(m3));  
  v=e(4)-e(3)-e(2)-e(1);
end

function m=ptranspose(m)
  scramble = [ 1 5 3 7 2 6 4 8 9 13 11 15 10 14 12 16];
  m=reshape(m(scramble),4,4);
end