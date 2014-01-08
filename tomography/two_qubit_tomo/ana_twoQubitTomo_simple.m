function val = ana_twoQubitTomo_simple(file, config)
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
if isempty(file)
    return
end
file=sort(file);

if ~exist('config','var')
    config=struct();
end
% set up default config stuff
config = def(config, 'opts','fullplot concurrence cal singleQ ');   % Random boolean options
config = def(config, 't1Time',{'before', 'before'}); % defaults to try to get t1 from before {left, right}
config = def(config, 'Ymult', []); % default Y-multiplier = 1; can set to -1 if not calibrating tomo. empty=auto
%config = def(config, 'mfit','fid norm sqn'); % Things to fit in mfit mode.
config = def(config, 'rng', [-inf inf]);  % Range
config = def(config, 'calFile', '');  
config = def(config, 'figind', 200);
figind = config.figind;
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
    [dataexp covariances scantime s, grps]=ana_TomoUnpack_2(file,config);  
    if isopt(config, 'bellfind') || isopt(config, 'Bellfind')
      [bestangles, beststates, bestpauli, bestfids_variances, bestfids, norms, norms_variances, statefunc, bpv] = fitBellStates(dataexp,covariances,stonly);
      val.bestangles = bestangles;
      val.bestates = beststates;
      val.bestpauli = bestpauli;
      val.bestfids_variance = bestfids_variances;
      val.bestfids = bestfids;
      val.norms = norms;
      val.norm_variances = norms_variances;
      val.statefunc = statefunc;
      val.bpv = bpv;
   end
    save(savefile);
end

xvals = 1:size(dataexp,2);
xlab = 'xval'; %FIXME
val.grps = grps;

% if 1
% figure(200+figind); 
%     figind=figind+1;
%     clf; hold on;
%     colors=[1 0 0; 0 0 1 ; .7 0 .7];
%     cind=[1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 ];
%     l=size(dataexp,1);
%     for i=1:l
%         area(xvals, bestpauli(:,i)+2*(l-i), 2*(l-i),'FaceColor',colors(cind(i),:));
%     end
%     set(gca,'YTick',2*(0:1:(l-1)));    
%     set(gca,'YTickLabel',bs(end:-1:1));
%     set(gca,'YLim',[-1 l*2-1]);
%     set(gcf,'Name','Pauli Plot');       
% end

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

if isopt(config,'singleQ')
  figure(200+figind); clf;  hold on;
  figind=figind+1;
  weight1 = zeros(1,size(dataexp,2));
  weight2 = weight1;
  weight3 = weight1;
  for j=1:size(dataexp,2)
      weight1(j) = sum(dataexp(1:3,j).^2);
      weight2(j) = sum(dataexp(4:6,j).^2);
      weight3(j) = sum(dataexp(7:end,j).^2);
  end
  c1  = [1 0 0]; c2 = [0 0 1];
  plot(xvals,weight1,'Color',c1);
  plot(xvals,weight2,'Color',c2);
  plot(xvals,weight3,'Color',.5*(c1+c2));
  legend({'Qubit 1','Qubit 2','Two-Qubit-Weight'});
  xlabel(xlab);
  ylabel('Norm');
  title('Single Qubit Norms');
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
    figs = 200+(config.figind:figind-1);
    ppt=guidata(pptplot);
    set(ppt.e_file,'String',file{1});
    set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
    set(ppt.e_title,'String',regexprep(regexprep(file{1},'(sm_)|(\.mat)',''),'.*\',''));
    %set(ppt.e_title,'String','ana_twoQubitTomo');
    if  isopt(config,'cal')
        set(ppt.e_body,'String',[sprintf('two qubit state tomography (corrected)\n')]);
    else
        set(ppt.e_body,'String',[sprintf('two qubit state tomography (not corrected)\n')]);
    end
end

val.file = file;
val.data =dataexp;
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