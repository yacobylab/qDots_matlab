function [figs pars sdata ]=ana_echo(file,config)
%function [figs pars sdata] = ana_echo(file,config) 
% config is a config struct describing what to do.  File can be blank
% ts is the total evolution times as a vector, ie. 0.05:0.05:2.  If ts
%   is nan(default), try to guess from group xvals.
% opts can be ramsey, noscale, fitdecay, phase, period, frq, amp, guessxval, even,
%  afitdecay -- fit decay only if amplitude is bigger than 5%
% odd, ramq, pdiff, nocenter, nocolor, ramt2, fitoffset (for frq decay)
% t1 for right is ~=0.05
% t1 for left is ~0
% dbz is a list of which channels have dbz data. defaults to any pulsegroup w/ dbz in the name.
% rng is points to fit. [1 10] means first 10 xvals.  [5 inf] means 5 onwared.
% grng is groups to fit.  [ .5 inf] means .5 onward.  uses group xvals
% channel is a list of channels to fit.
% nodbz -- skip dbz fit
% offset - offset between lines in figure 100

figs=[];
pars = [];
plotnum=1;
fitdescr='';
persistent lastname;
persistent sdata_cache;
persistent fileinfo;
if ~exist('file','var') || isempty(file)
  if ~isstr(lastname)
      lastname='';
  end
  file=uigetfile('sm*.mat','ana_echo',lastname);
end

if isempty(file)
    return;
end
cache=0;
if strcmp(lastname,file) && ~isempty(sdata_cache) && ~isempty(fileinfo)
    st=dir(file);
    if st.bytes == fileinfo.bytes && st.datenum == fileinfo.datenum
        cache=1;
    end
    fileinfo=st;
end
lastname=file;

if ~exist('config','var')
    config=struct();
end

% Load the data file
if cache
  if ~isopt(config,'quiet') % bug, config not defined yet. 
    fprintf('Using cached copy of file %s\n',file);
  end
  sdata=sdata_cache;
else
  sdata=load(file);  % Load the scan, data.. This lets us auto-generate some options.
  fileinfo=dir(file);
  sdata_cache=sdata;
end
  data=sdata.data;
  scan=sdata.scan;

scantime=getscantime(scan,data);
pars.scantime = scantime;
% Parse options
config = def(config,'opts','');   % Random boolean options
config = def(config,'rng',[]);    % Range of points to fit, ie [10 inf];
config = def(config,'grng',[]);   % Group range points to fit, ie [1.1 1.5];
config = def(config,'channel',1); % Range of channels to fit. fixme to be smarter
config = def(config,'t1',nan);   % Default ratio of t1 to tmeas.  fixme to be smarter.
config = def(config,'frames',[]); % Default frames to fit.
%config = def(config,'xlabel','T (\mus)'); % Default xlabel for inter-group series (put into
%the ramsey section below) 
config = def(config,'dxlabel','T (ns)');  % Default xlabel for intra-group series
config = def(config,'fitopts','interp');  % Interpolate xbals for fits\
config = def(config,'fb',100);    % Base number of figures to output.
config = def(config,'spsize',[2 2]); % Number of subplots for 'minor' figures.
% cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]'))
config = def(config,'dbz', find(cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]')))); % Is there a dBz group?
config = def(config,'side',[]); % Side of dot examined.
config = def(config, 'ts', nan); % Group xvals.
config = def(config,'acut',nan);  % Cutoff amplitude

if isopt(config,'ramsey') || isopt(config,'Ramsey')
    if isopt(config, 'singlegroup')
        config.opts = [config.opts 'fitdecay nocenter gauss'];
    else
        config.opts=[config.opts 'amp per frq gauss fitdecay nocenter ramq ramt2 epsRMS'];        
    end
    config = def(config,'xlabel','eps (mV)');  % Default xlabel for intra-group series
else
    config = def(config,'xlabel','T (\mus)'); % Default xlabel for inter-group series
end
       
if isopt(config,'echo')
   config.opts=[config.opts 'freq amp per'];
end

if isempty(config.frames)
    config.frames=1:size(data{config.channel},1);
end
if isempty(config.side)
    switch scan.loops(1).getchan{1}
        case 'DAQ2'   
          config.side='right';
        case 'DAQ1'
          config.side='left';
        otherwise
          error('Unable to determine side');
    end
end
   
if isnan(config.t1)
    [t1t config.t1] = att1(config.side,scantime,'before');
end

notdbz=setdiff(1:length(scan.data.pulsegroups),config.dbz);

%======================
% Find the dt's for individual scan lines.  
for i=1:length(scan.data.pulsegroups)
    try
    xvt=plsinfo('xval', scan.data.pulsegroups(i).name, [],scantime);
    dxvt=diff(xvt,[],2) ~= 0;
    [mc ind]=max(sum(dxvt,2));    
    dt(i,:)=xvt(ind,:);
    if ismember(i,notdbz)
      xv(:,i)=xvt(:);
    end
    catch
        warning('plsinfo is corrupt. params may be meaningless');
        dt(i,:) = 1:size(data{config.channel},3);
    end
end
config = def(config, 'dt', dt); % dts
dt = config.dt; %hack to keep backward compatibility
if any(size(dt)==1)
   dt = repmat(dt,size(data{1},2),1); 
end
if isnan(config.ts)
% Guess the group xval from the params
if isopt(config,'guessxval') 
    pulseparams=[];    
    for i = notdbz
        p=plsinfo('params', scan.data.pulsegroups(i).name, [],scantime);
        if isempty(p)
            fprintf('Warining; no parameters on group %s\n',scan.data.pulsegroups(i).name);
            fprintf('Probable logging error.  Throwing out data.\n');
            p=repmat(nan,size(pulseparams,2),1);
        end
        pulseparams(:,i)=p;
    end
    % Find the varpar that changes
    [r c]=find(diff(pulseparams(:,notdbz),[],2) ~= 0);
    switch mode(c)
        case 2
            xlab = 'T $(\mu{}s)$';
        case 3
            xlab = '$\epsilon$ (mV)';
        otherwise
            xlab = '?';
    end
    if isnan(mode(c))
        ts=1:length(notdbz);
    else
      ts=pulseparams(mode(c),:)';
    end
else
    ts=[];
    [r c] = find(diff(xv,[],2) ~= 0);
    if ~isempty(r)        
      ts = xv(mode(r),:)';
    else
      fprintf('No xval changes from group to group.  Try guessxval\n');           
      ts = xv(1,:);
    end
end
else
    ts=config.ts;
    if (length(ts) ~= length(scan.data.pulsegroups))
        error('The length of TS (%d) must be the same as the number of groups (%d)\n',length(ts),length(scan.data.pulsegroups));
    end
end

ts=ts(:);
if isopt(config,'even')
  notdbz=notdbz(2:2:end);
end
if isopt(config,'odd')
  notdbz=notdbz(1:2:end);
end
% Scale the data
if ~isopt(config,'noscale')
    %the next line is a hack. will break when ana_echo can do more than one
    %channel at once. 
  t1=ones(config.channel,1).*config.t1;
  %data_all=anaHistScale(scan,data,config.t1);
  data_all=anaHistScale(scan,data,t1);
  yl='P(T)';
  if isfield(config,'offset')
    offset=config.offset;
  else
    offset=1/3;
  end
else
  data_all=data;
  yl='V_{rf} (mV)';
  offset=4e-4;
end

if isopt(config, 'offsetoff')
   offset = 0; 
end

%pars.pulseparams = pulseparams;
omega_dbz=2*pi/32; % Assume this if there is no dbz reference.

for i=1:length(config.channel)
data=data_all{config.channel(i)};
fb=config.fb+100*(i-1);
 
  %=================
  % Plot the main data
  o=0;
  if ~isopt(config,'noplot')
    figure(fb); figs=[figs gcf];
    clf;  
  end
  colors='rgbcmk';
  color=@(x) colors(mod(x,end)+1);
  params=[];
  for j=1:length(notdbz) % Fit all the rest of the data.
     ind=notdbz(j);
     if isopt(config, 'singlegroup')
        rdata = squeeze(nanmean(data(config.frames,:))); 
     else
        rdata=squeeze(nanmean(data(config.frames,ind,:),1))'+o;
     end
     o=o+offset;
     if ~isopt(config, 'noplot')
         if isopt(config,'lines')
           plot(dt(ind,:),rdata,[color(j) '.-']);
         else
           plot(dt(ind,:),rdata,[color(j) '.']);
         end
     end
     if ~isopt(config,'nofit')
       if isopt(config,'noplot')
         [fp,ff]=fitosc(dt(ind,:),rdata,[config.opts],config.rng,[color(j) '-']);
       else
         [fp,ff]=fitosc(dt(ind,:),rdata,['plot' config.opts],config.rng,[color(j) '-']);
       end           
       params(j,:)=fp;     
     end
     if ~isopt(config,'noplot')
       ylabel(yl);
       xlabel(config.dxlabel);
       hold on;
     end
  end
  if isopt(config,'nofit')
      return;
  end
  pars.params=params;
  pars.ts=ts;
  pars.dt=dt;
  if length(params) > 3
      if size(params,1) > 10
        jbar=mean(params(3:end-3,4));
      else
        jbar=mean(params(:,4));
      end
  else
      jbar=nan;
  end
  title(sprintf('|J|=%g Mhz',1e3*jbar/(2*pi)));
  fitdescr = [ fitdescr sprintf('|J|=%g Mhz\n',1e3*jbar/(2*pi)) ];
  
  if ~isopt(config,'nodbz')
      for i=1:length(config.dbz)
          dbzdata=squeeze(nanmean(data(config.frames,config.dbz(i),:),1))';
          [plotnum figs] = nextfig(config,plotnum,fb,figs);
          plot(dt(1,:),dbzdata,'b.');
          xlabel('T (ns)');
          ylabel(yl);
          
          if ~isopt(config,'nofitdbz')
              [fp,ff]=fitosc(dt(1,:),dbzdata,['fitdecay nocenter plot ' config.fitopts],[]);
              hold on;
              str=sprintf('T_2^*=%.3g ns, V=%.3f, T=%.3f, phi=%f',1./fp(6),2*sqrt(fp(2)^2+fp(3)^2),2*pi/fp(4),atan2(fp(3),fp(2))-pi);
              title(str);
              pars.dbzt2=1./fp(6);
              fitdescr = [ fitdescr sprintf('dBz_%d: ',config.dbz(i)) str sprintf('\n') ];
              omega_dbz = fp(4);
              pars.omega_dbz=omega_dbz;
          else
              title('dBz reference');
          end
      end
  end
  
  
  ts=ts(notdbz);
  % Plot various handy things.
  % amplitude vs. xval
  if isopt(config,'amp')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);
    %ampfunc=@(x) 2*((abs(x(:,4)) > 0.01) .* (abs(x(:, 6)) < .2).*(sqrt(x(:, 2).^2 + x(:, 3).^2)));  
    ampfunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));  
    a=ampfunc(params);
    plot(ts,a,'b.');   
    if ~isnan(config.acut)       
       cut=config.acut*median(a(1:3));
       mask=(a > cut);
    else
       mask=~isnan(a);
    end
    if isopt(config, 'rmoutlier')
        mask = mask & (a <1.2);
    end
    plot(ts(mask),ampfunc(params(mask,:)),'b.');   
    if isopt(config,'linfit')
      ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask);
      fpp = fitwrap('plfit',ts(ind)', ampfunc(params(ind,:))', [-1 0], @(p,x) p(1)*x + p(2));
      str = [str sprintf('Amp %g t + %g',fpp(1),fpp(2))];
      fp=[];
      fp(1) = fpp(2); fp(2) = (fpp(2)+fpp(1)*mean(ts(ind)))/fpp(1);
    elseif isopt(config,'logfit')
      ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask(:)');
      fpp = fitwrap('plfit',ts(ind)', log(ampfunc(params(ind,:)))', [-1 0], @(p,x) p(1)*x + p(2));
      str = [str sprintf('Amp exp^(%g t + %g)',fpp(1),fpp(2))];
      fp=[];
      fp(1) = fpp(2); fp(2) = 1/fpp(1);
    else
      [fp,junk,fstr]=fitdecay(ts(mask)',ampfunc(params(mask,:))',['plot' config.opts],config.grng);
      str=['Amp' fstr];
      str= [str sprintf('Amp=%.3g T_2^{echo}=%.3g Q=%.1f T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/jbar),2*pi/mean(jbar))];
    end
    title(str);
    fitdescr = [ fitdescr 'Amp: ' str sprintf('\n') ];
    ylabel(yl);
    xlabel(config.xlabel);
    pars.amp = fp(1);
    if ~isempty(config.grng)
      ind = find(ts > config.grng(1) & ts < config.grng(2));
    else
      ind=1:length(ts);
    end 
    pars.maxamp=max(ampfunc(params(ind,:)));
    pars.T2 = fp(2);
    pars.Q = fp(2)/(1e-3*2*pi/jbar);
    pars.Jbar = jbar;
    pars.T= 2*pi/mean(jbar);
    pars.afp = fp;
  end
  
  % Plot various handy things.
  % amplitude vs. xval
  if isopt(config,'per')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);
    perfunc=@(params) 2*pi./params(:,4);
    plot(ts,perfunc(params),'b.');
    title('Period');
    ylabel('T (ns)');
    xlabel(config.xlabel);
  end
  
  if isopt(config,'frq')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);    
    if ~isopt(config,'noref')
      freqfunc=@(params) sqrt( abs(params(:,4)./(2*pi)).^2-((omega_dbz/((2*pi)))^2)) .* sign(params(:,4) - omega_dbz);   
    else
      freqfunc=@(params) params(:,4)./(2*pi);
    end
    plot(ts,1e3*freqfunc(params),'b.');
    %fp=fitdecay(ts',ampfunc(params)',['plot' opts]);
    %title(sprintf('Amp=%.3g T_2^{echo}=%.3g Q=%.1f T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/jbar),2*pi/mean(jbar)));
    ylabel('J (MHz)');
    xlabel(config.xlabel);
    if isopt(config,'fitoffset')
      [decp decf] = fitdecay(ts',1e3*freqfunc(params)','plot fitoffset',config.grng);
    else
      [decp decf] = fitdecay(ts',1e3*freqfunc(params)','plot',config.grng);
    end
    str=sprintf('decay const = %.2d', decp(2));    
    title(str);
    fitdescr = [ fitdescr 'Freq: ' str sprintf('\n') ];
    % 2nd form is more useful for quickly guestimating J's.
    %fitdescr = [ fitdescr sprintf('J(eps) = %.3g * exp(-eps/%.3g) + %.3g Mhz',decp(1:3)) ];
    fitdescr = [ fitdescr sprintf('J(eps) = 100 * exp(-(eps-%.4g)/%.3g) + %.3g Mhz',log(decp(1)/100)*decp(2),decp(2),decp(3))] ;
    pars.freqfunc=@(eps) 100*exp(-(eps-log(decp(1)/100)*decp(2))/decp(2))+decp(3);
    pars.decp = decp;
    
    pars.freq = [ts';1e3*freqfunc(params)'];
  end
    if isopt(config,'ramt2')
       [plotnum figs] = nextfig(config,plotnum,fb,figs);
       plot(ts,abs(1./params(:,6)));
       xlabel(config.xlabel);
       ylabel('T_2^*, ns');
    end
    
    if isopt(config,'epsRMS')|| isopt(config,'epsrms')
       [plotnum figs] = nextfig(config,plotnum,fb,figs);
       %dJdE is in MHz/mV
       djde = (1/decp(2))*(1e3*freqfunc(params)-decp(3)); %Offset does not contribute to slope
       plot(1./djde(ind),abs(1./params(ind,6)),'.'); hold on; % T2* here is in ns
       noiseslope = (1./djde(ind))'/abs(1./(params(ind,6)))'; % this is matlab shorthand for least squares slope
       plot([0 1/(djde(ind(end)))],[0 1/(djde(ind(end)))]./noiseslope,'g');
%       erms = noiseslope/(2*pi*sqrt(2)); % No mikey; see our paper pg. 4
       erms=sqrt(2)*noiseslope/(2*pi);     
% dJ/dE is in MHz/mV=GHz/V, t2* is in ns, so erms in in v.
       xlabel('(dJ/deps)^{-1} (mV/MHz)');
       ylabel ('T2^* (ns)');
       title(sprintf('RMS Voltage noise =%g (\\muV)',erms*1e6));
       fitdescr = [ fitdescr sprintf('RMS Noise: %g uV\n',erms*1e6)];       
       [plotnum figs] = nextfig(config,plotnum,fb,figs);
       plot(1e3*freqfunc(params),djde);
       xlabel('J');
       ylabel('dJ/d\epsilon (MHz/mV)');
    end
  if isopt(config,'echocenter')
      [plotnum figs] = nextfig(config,plotnum,fb,figs);
      plot(ts,params(:,5));
      xlabel(config.xlabel)
      ylabel('Echo center (ns)');
  end
  if isopt(config,'echophase')
      [plotnum figs] = nextfig(config,plotnum,fb,figs);
      plot(ts,unwrap(atan2(params(:,3),params(:,2))));
      xlabel(config.xlabel)
      ylabel('Echo Phase (radians)');
  end      
  if ~isopt(config,'nocolor')
     [plotnum figs] = nextfig(config,plotnum,fb,figs);
     rdata=reshape(permute(data,[1 3 2]),size(data,1),size(data,2)*size(data,3));     
     imagesc(rdata(config.frames,:));
  end
  if isopt(config,'mean')
     [plotnum figs] = nextfig(config,plotnum,fb,figs);
     rdata=nanmean(data(config.frames,:,:),1);     
     imagesc(squeeze(rdata));
  end
  
  if isopt(config,'echonoise') && isfield(config,'djde') && ~isempty(config.djde)
      %this only works with power and exponential decay. not 'both'
      djde = config.djde*1e9; %MHz/mV into Hz/V
      if isopt(config,'power')
          beta = fp(4)-1;
          gval=gamma(-1-beta)*sin(pi*beta/2);
      else
          beta = 0;
          gval=pi/2;  % Limit of gval above as beta->0
      end
      pars.Se = (2*pi/abs(2^-beta*(-2+2^beta)*(1e-6*pars.T2)^(1+beta)*(2*pi*djde)^2*gval));
      % In limit as beta->0, this is 1/(4*pi^2)
      pars.Se = pars.Se / ( (2*pi)^beta); % 1/f, not 1/omega
      pars.Seps = @(f) pars.Se/f^beta;
      fitdescr = [ fitdescr sprintf('Noise@1Mhz: %g nV (beta=%g)\n',sqrt(pars.Seps(1e6))*1e9,beta)];      
  end
  
end
%pars = params;
figs = unique(figs);

  if ~isopt(config,'noppt')
     prettyfile=regexprep(file,'(sm_)|(\.mat)','');
     indentdescr= regexprep(fitdescr,'^(.)','\t$1','lineanchors');
     ppt=guidata(pptplot);     
     set(ppt.e_file,'String',file);     
     set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);     
     set(ppt.e_title,'String',prettyfile);
     set(ppt.e_body,'String',fitdescr);     
     clipboard('copy',sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
     set(ppt.exported,'Value',0);
  end
return;

%   if ~isempty(strfind(opts,'phase'))
%       pl.col = name;
%       pl.x = 7;
%       pl.y='@(x)unwrap(atan2(x(:,2),x(:,3)))';      
%       pl.params=ts';
%       f3=makefits(col, pl, ds(1), length(ds));
%       dvdisplay(col,f3,'flag9');
%       dvdisplay(col,f3,'title','String','Ramsey Phase');
%       dvdisplay(col,f3,'xlabel','String',xlab,'interpreter','latex');
%       dvdisplay(col,f3,'ylabel','String','Phase (rad)');
%       dsa=union(dsa,f3);
%   end
%     
%   if ~isempty(strfind(opts,'pdiff'))
%       pl.col = name;
%       pl.x = '@(x) x(1:2:end,7)';
%       pl.y='@(x)unwrap(atan2(x(1:2:end,2),x(1:2:end,3))-atan2(x(2:2:end,2),x(2:2:end,3)))';
%       pl.params=ts';
%       f3=makefits(col, pl, ds(1), length(ds));
%       dvdisplay(col,f3,'flag9');
%       dvdisplay(col,f3,'title','String','Ramsey Phase Change');
%       dvdisplay(col,f3,'xlabel','String',xlab,'interpreter','latex');
%       dvdisplay(col,f3,'ylabel','String','Phase (rad)');
%       dsa=union(dsa,f3);
%       
%       dvfit(col,f3,'plint plfit ause',[1 1], '@(p,x) p(1)*x+p(2)',[1 1 ]);
%       fp=dvfit(col,f3,'getp');
%       dvdisplay(col,f3,'flag9');  % Seems to mean 'title'; talk to Hendrik
%       dvdisplay(col,f3,'title',...
%         'String',sprintf('$d\\phi/dt = %.3f rad/\\mu{}sec$',fp(1)),...
%         'Interpreter','latex');
%     
%     
%       pl.col = name;
%       pl.x = '@(x) x(1:2:end,7)';
%       pl.y='@(x)unwrap(atan2(x(1:2:end,2),x(1:2:end,3))-atan2(x(2:2:end,2),x(2:2:end,3)))';
%       pl.params=ts';
%       f3=makefits(col, pl, ds(1), floor(length(ds)/4));
%       dvdisplay(col,f3,'flag9');
%       dvdisplay(col,f3,'title','String','Early Ramsey Phase Change');
%       dvdisplay(col,f3,'xlabel','String',xlab,'interpreter','latex');
%       dvdisplay(col,f3,'ylabel','String','Phase (rad)');
%       dsa=union(dsa,f3);
%       
%       dvfit(col,f3,'plint plfit ause',[1 1], '@(p,x) p(1)*x+p(2)',[1 1 ]);
%       fp=dvfit(col,f3,'getp');
%       dvdisplay(col,f3,'flag9');  % Seems to mean 'title'; talk to Hendrik
%       dvdisplay(col,f3,'title',...
%         'String',sprintf('$d\\phi/dt = %.3f rad/\\mu{}sec$',fp(1)),...
%         'Interpreter','latex');
% 
%   end
% 
%  
%      
% 
%   if ~isempty(strfind(opts,'ramt2'))
%       pl.col = name;
%       pl.x=7;
%       pl.y = '@(x) 1./x(:,6)';
%       pl.params=ts';
%       f3=makefits(col, pl, ds(1), length(ds));
%       dsa=union(dsa,f3);
%       dvdisplay(col,f3,'flag9');  % Seems to mean 'title'; talk to Hendrik
%       dvdisplay(col,f3,'title',...
%           'String',sprintf('Ramsey T_2^*',fp(1),fp(2)),...
%           'Interpreter','tex');
%       dvdisplay(col,f3,'ylabel','String','T_2^*','interpreter','latex');
%       dvdisplay(col,f3,'xlabel','String',xlab,'interpreter','latex');
%   end
% 
% end
%    dvplot(col,dsa);    
% return;
% 
% % Make sure there are no scoping issues
% function fitdecay(col,f)  
%   dvfit(col,f,'plint plfit ause',[1 1 0], '@(p,x) p(1)*exp(-x/p(2))+p(3)',[1 1 0]);
% return;
% 

function [fp,ff]=fitosc(x,y,opts,rng,style)
   % initialization function
   fig=gcf;
   fifn.fn = @fioscill;
   fifn.args = {1};
 
   cosfn = '@(y, x)y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)';
   cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
   cosfn3 = '@(y,x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*(y(5)*x)'; %linear decay
   cosfn4 = '@(y, x)y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)';
   cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
      
   if exist('rng','var') && ~isempty(rng)
       pts=x>rng(1) & x <rng(2);
       x=x(pts);
       y=y(pts);
   end

   if isempty(strfind(opts,'fitdecay')) || (~isempty(strfind(opts,'afitdecay')) && std(y) < 2e-2)     
     fifn.args={2}; % No decay
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     ig = [fp(1), fp(2)*cos(fp(3)), fp(2)*sin(fp(3)), fp(4:6)];
     fp=fitwrap('fine',x,y,ig, cosfn, [1 1 1 1 0 0]);
     ff=str2func(cosfn);         
   elseif isempty(strfind(opts,'nocenter'))
     fifn.args={2}; % Decay and center
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     ig = [fp(1), fp(2)*cos(fp(3)), fp(2)*sin(fp(3)), fp(4:6)];
     fp=fitwrap('fine',x,y,fifn,cosfn2, [1 1 1 1 0 0]);
     fp=fitwrap('fine',x,y,fp, cosfn2, [1 1 1 1 1 1]);
     ff=str2func(cosfn2);
   else      % Decay but no center
     fifn.args={2};
     fp=fitwrap('fine',x,y,fifn, cosfn5, [1 0 1 1 0 0]);
     fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
     fp=fitwrap('fine',x,y,fp, cosfn2, [1 1 1 1 0 1]);
     ff=str2func(cosfn2);
   end
   if ~isempty(strfind(opts,'plot'))
       figure(fig);
       hold on;
       if ~isempty(strfind(opts,'interp'))
           x=linspace(x(1,1),x(1,end),512);
       end
       if exist('style','var') && ~isempty('style')
         plot(x,ff(fp,x),style);
       else
         plot(x,ff(fp,x),'r-');
       end
   end
return;

function [fp,ff,fitstring]=fitdecay(x,y,opts,rng,style)
   % initialization function
   fig=gcf;

   if isempty(strfind(opts,'fitoffset'))
      mask=[1 1 0];
   else
      mask=[1 1 1];
   end
      
   if exist('rng','var') && ~isempty(rng)
       pts=x>rng(1) & x <rng(2);
       x=x(pts);
       y=y(pts);
   end

   if ~isempty(strfind(opts,'gauss'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^2)+p(3);
     init=[1 max(x)/3 min(y)];
     init(1)=range(y)/(ff(init,min(x)-init(3)));     
     fmt='Decay: %.3f exp(-(t/%3f)^2)+%g\n';
     fperm=[1 2 3];
   elseif ~isempty(strfind(opts,'both'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^2-x/p(4))+p(3);
     init=[1 max(x)/3 min(y) max(x)/10];
     init(1)=range(y)/(ff(init,min(x))-init(3));
     mask(4)=1;
     fperm=[1 2 3];
     fmt='Decay: %.3f exp(-t/%3f -(t/%3f)^2)+%g\n';
   elseif ~isempty(strfind(opts,'power'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^p(4))+p(3);       
     init=[1 max(x)/3 min(y) 1];
     mask(4)=1;
     init(1)=range(y)/(ff(init,min(x))-init(3));   
     fperm=[1 2 4 3];
     fmt='Decay: %.3f exp(-(t/%3f)^{%f})+%g\n';   
   else       
     ff=@(p,x) p(1)*exp(-x/p(2))+p(3);       
     init=[1 max(x)/3 min(y)];
     init(1)=range(y)/(ff(init,min(x))-init(3));
     fmt='Decay: %.3f exp(-t/%3f)+%g\n';
     fperm=[1 2 3];
   end

   if isempty(strfind(opts,'fitoffset'))      
      init(3)=0;
   end

   fp = fitwrap('plinit plfit',x,y,init,ff,mask);
   fitstring=sprintf(fmt,fp(fperm));
   if ~isempty(strfind(opts,'plot'))
       figure(fig);
       hold on;
       if ~isempty(strfind(opts,'interp'))
           x=linspace(x(1,1),x(1,end),512);
       end
       if exist('style','var') && ~isempty('style')
         plot(x,ff(fp,x),style);
       else
         plot(x,ff(fp,x),'r-');
       end
   end
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

function [plotnum figs] = nextfig(config, plotnum, fb, figs)   
   figure(1+fb+floor((plotnum-1)/prod(config.spsize))); 
   if(mod(plotnum-1,prod(config.spsize)) == 0)
       clf;
   end   
   subplot(config.spsize(1),config.spsize(2), mod(plotnum-1,prod(config.spsize))+1); 
   figs=unique([figs gcf]); 
   plotnum=plotnum+1;
return